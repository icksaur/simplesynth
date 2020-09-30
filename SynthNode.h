#pragma once

#include "synth.h"
#include <vector>

const StereoSample emptySample{};

StereoSample operator+(StereoSample a, StereoSample b) {
    return StereoSample{ a.l + b.l, a.r + b.r };
}

StereoSample operator*(StereoSample a, float b) {
    return StereoSample{ a.l * b + a.r * b };
}

struct SynthParameter {
    StereoSample setting;
    StereoSample accumulator;
    const char* name;
};

struct SynthRoute {
    class SynthNode* target;
    float multiplier;
    int parameterIndex;
};

#define PARAM_COUNT 6

struct SynthNode {
    std::vector<SynthRoute> routes;
    std::vector<SynthParameter> inputs;
    SynthNode() = delete;
    SynthNode(int inputs) : inputs(inputs) {

    }
    void set(int inputIndex, float settingValue) {
        inputs[inputIndex].setting = StereoSample{ settingValue };
    }
    void addRoute(const SynthRoute route) {
        routes.push_back(route);
    }
    void accumulate(int index, StereoSample value) {
        inputs[index].accumulator = inputs[index].accumulator + value;
    }
    float stepFloatInput(int index) {
        float accumulated = inputs[index].accumulator.l;
        inputs[index].accumulator = emptySample;
        return inputs[index].setting.l + accumulated;
    }
    StereoSample stepSampleInput(int index) {
        StereoSample accumulated = inputs[index].accumulator;
        inputs[index].accumulator = emptySample;
        return inputs[index].setting + accumulated;
    }
    virtual StereoSample stepInternal() = 0;
    void step() {
        StereoSample result = stepInternal();
        for (SynthRoute& route : routes) {
            route.target->accumulate(route.parameterIndex, result * route.multiplier);
        }
    }

};

struct LFONode : public SynthNode {
    LFO lfo;
    LFONode() : SynthNode(2), lfo{} {
        set(0, 1.0f);
        set(1, 0.5f);
    }
    StereoSample internalStep() {
        lfo.amplitude = stepFloatInput(0);
        lfo.period = stepFloatInput(1);
        return StereoSample{ lfo.generate() };
    }
};

struct OscillatorNode : public SynthNode {
    Note noteGen;
    Oscillator oscillator;
    OscillatorNode() : SynthNode(6), noteGen{}, oscillator{} {
        set(0, 0.001f);
        set(1, 4.0f);
        set(2, 1.0f);
    }
    StereoSample internalStep() {
        noteGen.note = stepFloatInput(0);
        noteGen.glide = stepFloatInput(1);
        oscillator.sine = stepFloatInput(2);
        oscillator.saw = stepFloatInput(3);
        oscillator.square = stepFloatInput(4);
        oscillator.noise = stepFloatInput(5);
        float note = noteGen.generate();
        float freq = noteFrequency(note);
        float sample = oscillator.generate(freq);
        return StereoSample{ sample };
    }
};

struct EnvelopeNode : public SynthNode {
    Envelope envelope;
    EnvelopeNode() : SynthNode(4), envelope{ 0.05f, 0, 1, 0.15f, 0, 0, 4 } {
        
    }
    StereoSample internalStep() {
        float gate = stepFloatInput(0);
        envelope.a = stepFloatInput(1);
        envelope.d = stepFloatInput(2);
        envelope.s = stepFloatInput(3);
        envelope.r = stepFloatInput(4);
        return StereoSample{ envelope.generate(gate) };
    }
};

struct FilterNode : public SynthNode {
    MoogFilter filter;
    FilterNode() : SynthNode(3), filter(0.5f, 0.0f) {
        inputs[1].setting = StereoSample{ 0.5f };
        inputs[2].setting = StereoSample{ 0.0f };
    }
    StereoSample internalStep() {
        filter.set(stepFloatInput(1), stepFloatInput(2));
        float sample = stepFloatInput(0);
        filter.process(sample);
        return StereoSample{ sample };
    }
};

struct VocalNode : public SynthNode {
    VocalFilter filter;
    VocalNode() : SynthNode(1) {}
    StereoSample internalStep() {
        float sample = stepFloatInput(0);
        filter.process(sample);
        return StereoSample{ sample };
    }
};

struct Reverb : public SynthNode {
    StereoReverb reverb;
    Reverb() : SynthNode(3) {
        inputs[1].setting = StereoSample{ 0.82f };
        inputs[2].setting = StereoSample{ 0.2f };
    }
    StereoSample internalStep() {
        StereoSample result;
        reverb.set(stepFloatInput(1), stepFloatInput(2));
        reverb.process(stepFloatInput(0), result);
    }
};

struct DelayNode : public SynthNode {
    Delay delay;
    DelayNode() : SynthNode(3), delay{} {
        set(1, 0.5f);
        set(2, 0.5f);
    }
    StereoSample internalStep() {
        float input = stepFloatInput(0);
        delay.delay = stepFloatInput(1);
        delay.decay = stepFloatInput(2);
        delay.process(input);
        return StereoSample{ input };
    }
};

struct StereoNode : public SynthNode {
    Fader fader;
    StereoNode() : SynthNode(2), fader{} {
        set(1, 0.0f);
    }
    StereoSample internalStep() {
        StereoSample sample = stepSampleInput(0);
        fader.process(stepFloatInput(1), sample);
        return sample;
    }
};

struct OutputNode : public SynthNode {
    StereoSample* buffer;
    size_t sampleIndex;
    OutputNode() : SynthNode(1), buffer(nullptr), sampleIndex(0) {}
    StereoSample internalStep() {
        buffer[sampleIndex++] = stepSampleInput(0);
        return emptySample;
    }
};
