#include "synth.h"
#include <cmath>
#include <cstring>
#include <algorithm>
#include <cstdlib>

// StereoSample operations
const StereoSample emptySample{};

StereoSample operator+(StereoSample a, StereoSample b) {
    return StereoSample{ a.l + b.l, a.r + b.r };
}

StereoSample operator*(StereoSample a, float b) {
    return StereoSample{ a.l * b, a.r * b };
}

// Formant vowel definitions
FormantVowel soprano_a = {800, 1150, 2900, 3900, 4950, 0, -6, -32, -20, -50, 80, 90, 120, 130, 140};
FormantVowel soprano_e = {350, 2000, 2800, 3600, 4950, 0, -20, -15, -40, -56, 60, 100, 120, 150, 200};
FormantVowel soprano_i = {270, 2140, 2950, 3900, 4950, 0, -12, -26, -26, -44, 60, 90, 100, 120, 120};
FormantVowel soprano_o = { 450, 800, 2830, 3800, 4950, 0, -11, -22, -22, -50, 70, 80, 100, 130, 135};
FormantVowel soprano_u = { 325, 700, 2700, 3800, 4950, 0, -16, -35, -40, -60, 50, 60, 170, 180, 200};
FormantVowel alto_a = { 800, 1150, 2800, 3500, 4950, 0, -4, -20, -36, -60, 80, 90, 120, 130, 140};
FormantVowel alto_e = { 400, 1600, 2700, 3300, 4950, 0, -24, -30, -35, -60, 60, 80, 120, 150, 200};
FormantVowel alto_i = { 350, 1700, 2700, 3700, 4950, 0, -20, -30, -36, -60, 50, 100, 120, 150, 200};
FormantVowel alto_o = { 450, 800, 2830, 3500, 4950, 0, -9, -16, -28, -55, 70, 80, 100, 130, 135};
FormantVowel alto_u = { 325, 700, 2530, 3500, 4950, 0, -12, -30, -40, -64, 50, 60, 170, 180, 200};
FormantVowel countertenor_a = { 660, 1120, 2750, 3000, 3350, 0, -6, -23, -24, -38, 80, 90, 120, 130, 140};
FormantVowel countertenor_e = { 440, 1800, 2700, 3000, 3300, 0, -14, -18, -20, -20, 70, 80, 100, 120, 120};
FormantVowel countertenor_i = { 270, 1850, 2900, 3350, 3590, 0, -24, -24, -36, -36, 40, 90, 100, 120, 120};
FormantVowel countertenor_o = { 430, 820, 2700, 3000, 3300, 0, -10, -26, -22, -34, 40, 80, 100, 120, 120};
FormantVowel countertenor_u = { 370, 630, 2750, 3000, 3400, 0, -20, -23, -30, -34, 40, 60, 100, 120, 120};
FormantVowel tenor_a = { 650, 1080, 2650, 2900, 3250, 0, -6, -7, -8, -22, 80, 90, 120, 130, 140};
FormantVowel tenor_e = { 400, 1700, 2600, 3200, 3580, 0, -14, -12, -14, -20, 70, 80, 100, 120, 120};
FormantVowel tenor_i = { 290, 1870, 2800, 3250, 3540, 0, -15, -18, -20, -30, 40, 90, 100, 120, 120};
FormantVowel tenor_o = { 400, 800, 2600, 2800, 3000, 0, -10, -12, -12, -26, 40, 80, 100, 120, 120};
FormantVowel tenor_u = { 350, 600, 2700, 2900, 3300, 0, -20, -17, -14, -26, 40, 60, 100, 120, 120};
FormantVowel bass_a = { 600, 1040, 2250, 2450, 2750, 0, -7, -9, -9, -20, 60, 70, 110, 120, 130};
FormantVowel bass_e = { 400, 1620, 2400, 2800, 3100, 0, -12, -9, -12, -18, 40, 80, 100, 120, 120};
FormantVowel bass_i = { 250, 1750, 2600, 3050, 3340, 0, -30, -16, -22, -28, 60, 90, 100, 120, 120};
FormantVowel bass_o = { 400, 750, 2400, 2600, 2900, 0, -11, -21, -20, -40, 40, 80, 100, 120, 120};
FormantVowel bass_u = { 350, 600, 2400, 2675, 2950, 0, -20, -32, -28, -36, 40, 80, 100, 120, 120};

float noteFrequency(float a) {
    return 27.5f * pow(2.0f, (float)a);
}

float noteFrequency(int a, int halfSteps) {
    return noteFrequency((float)a + (float)halfSteps / 12.0f);
}

void Envelope::trigger() {
    _gate = 1;
    state = 0; 
    t = 0;
}

float Envelope::generate(int gate) {
    float amp1;
    if (_gate && !gate) {
        if (state < 3) {
            state = 3; t = 0;
        }
        _gate = gate;
    }
    switch(state) {
    case 0:
        if(a && t < 1) {
            amp1 = t;
            t += 1 / a / (float)frequency;
            break;
        }
        t = 0;
        state++;
    case 1:
        if (d && state < 2) {
            amp1 = (1 - t * (1 - s));
            t += 1 / d / (float)frequency;
            if (t >= 1) {
                state++;
                t = 0;
            }
            break;
        }
        t = 0;
        state++;
    case 2:
        amp1 = s;
        break;
    case 3:
        if (r && t < 1) {
            amp1 = s * (1-t);
            t += 1 / r / frequency;
            break;
        }
        t = 0;
        state = 4;
    case 4:
        amp1 = 0;
    }
    float ampDif = amp1 - amp;
    if (fabs(ampDif) > envelopePopThreshold) {
        amp1 = amp + copysign(envelopePopThreshold, ampDif);
    }
    amp = amp1;
    return amp1;
}

float Note::generate() {
    float n = (note * t + note0 * (1 - t));
    t += 1 / (float)frequency / glide;
    if (t >= 1) {
        t = 1;
        note0 = note;
    }
    return n;
}

float Oscillator::generate(float f) {
    float h = sin(2.0f * (float)M_PI * t) * sine;
    h += (t / 0.5f - 1.0f) * saw;
    h += (t < 0.5 ? -1 : 1) * square;
    h += ((float)rand() / (float)RAND_MAX * 2.0f - 1.0f) * noise;
    t += f / frequency;
    wrap(t);
    return h;
}

float LFO::generate() {
    float val = (float)sin(2.0 * (float)M_PI * t) * amplitude;
    t += 1 / (float)frequency / period;
    wrap(t);
    return val;
}

StagedFilter::StagedFilter() {
    k = 0.5f;
    memset(y, 0, sizeof(y));
}

void StagedFilter::process(float & sample) {
    float k = fmin(this->k, 0.95f); // cache local K less than 0.5f
    float x = sample;
    sample = pow(1 - k, 4) * x;
    sample += (4.0f * k) * y[0];
    sample += (6.0f * k * k) * -y[1];
    sample += (4.0f * pow(k, 3)) * y[2];
    sample += pow(k, 4) * -y[3];
    
    y[3] = y[2];
    y[2] = y[1];
    y[1] = y[0];
    y[0] = sample;
}

BandpassFilter::BandpassFilter() {
    f = 0.2f;
    r = 0.95f;
    memset(x, 0, sizeof(x));
    memset(y, 0, sizeof(y));
}

void BandpassFilter::process(float & sample) {
    float x0 = sample;
    float b0 = (1.0f - r * r) / 2.0f;
    sample = sample * b0;
    //sample += x[1] * -b0; // improve response - better without this for vocal filter
    sample -= y[0] * -2.0f * r * cos(2.0f * (float)M_PI * f);
    sample -= y[1] * r * r;
    //improvement: multiply this by some -dB gain to balance formants
    
    y[1] = y[0];
    y[0] = sample;
    x[1] = x[0];
    x[0] = x0;
}

MoogFilter::MoogFilter(float frequency, float resonance) {
    b0 = b1 = b2 = b3 = b4 = t1 = t2 = 0.0f;
    type = low;
    set(frequency, resonance);
}

void MoogFilter::panic() {
    b0 = b1 = b2 = b3 = b4 = t1 = t2 = 0.0f;
}

void MoogFilter::set(float frequency, float resonance) {
    q = 1.0f - frequency;
    p = frequency + 0.8f * frequency * q;
    f = p + p - 1.0f;
    q = resonance * (1.0f + 0.5f * q * (1.0f - q + 5.6f * q * q));
}

void MoogFilter::process(float & sample) {
    float & in = sample;
    in -= q * b4;                          //feedback
    t1 = b1;  b1 = (in + b0) * p - b1 * f;
    t2 = b2;  b2 = (b1 + t1) * p - b2 * f;
    t1 = b3;  b3 = (b2 + t2) * p - b3 * f;
    b4 = (b3 + t1) * p - b4 * f;
    b4 = b4 - b4 * b4 * b4 * 0.166667f;    //clipping
    b0 = in;

    switch (type) {
    case low: sample = b4; break;
    case high: sample = in - b4; break;
    case band: sample = 3.0f * (b3 - b4); break;
    }
}

float MoogFreq(float f) {
    return exp(2.0f * (float)M_PI * f);
}

float peak(int hz) {
    return (float)hz / (float)frequency;
}

float bandwidth(int hz) {
    return 1.0f - (float)hz / (float)frequency;
}

float gain(int db) {
    return ((float)100 + (float)db) / 100.0f;
}

VocalFilter::VocalFilter() {
    formants = 5;
    setVowelFormants(&soprano_a);
}

void VocalFilter::setVowelFormants(FormantVowel * vowel) {
    for (int i = 0; i < 5; i++) {
        filters[i].f = peak(vowel->hz[i]);
        filters[i].r = bandwidth(vowel->bw[i]);
    }
}

void VocalFilter::process(float & sample) {
    float x0 = sample, x = sample;
    sample = 0.0f;
    for (int i = 0; i < formants; i++) {
        filters[i].process(x);
        sample += x;
        x = x0;
    }
}

void AllPassDelay::process(float & sample) {
    float r1 = buffer[read] + sample * k; // peek in delay
    sample = r1 * -k + buffer[read];
    buffer[(read + delay) % frequency] = r1 * 0.5f;
    ++read %= frequency;
}

void AllPassFilter::process(float & sample) {
    float r1 = sample + r * k;
    sample = r1 * -k + r;
    r = r1;
}

LowpassFeedbackCombFilter::LowpassFeedbackCombFilter(float room, float damp, int delay) {
    if (delay >= combSize) throw std::runtime_error("delay is too long");
    this->damp = damp;
    this->room = room;
    this->n = delay;
    buffer = new float[combSize]; // guess this could be n
    memset(buffer, 0, sizeof(float) * combSize);
    read = combSize - n;
    onepole = 0.0f;
}

LowpassFeedbackCombFilter::~LowpassFeedbackCombFilter() {
    delete[] buffer;
}

float LowpassFeedbackCombFilter::process(float x) {
    float y = buffer[read]; // delay line
    onepole = (1.0f - damp) / (1.0f - damp * onepole);
    int write = (read + n) % combSize;
    buffer[write] = x + room * onepole * y;
    read = (read + 1) % combSize;
    return y;
}

MonoReverb::MonoReverb() {
    combs[0] = new LBCF(defaultRoom, .2f, 1617);
    combs[1] = new LBCF(defaultRoom, .2f, 1557);
    combs[2] = new LBCF(defaultRoom, .2f, 1491);
    combs[3] = new LBCF(defaultRoom, .2f, 1277);
    combs[4] = new LBCF(defaultRoom, .2f, 1422);
    combs[5] = new LBCF(defaultRoom, .2f, 1356);
    combs[6] = new LBCF(defaultRoom, .2f, 1188);
    combs[7] = new LBCF(defaultRoom, .2f, 1116);
    allpass[1] = new AllPassDelay{allPassGain, 556};
    allpass[2] = new AllPassDelay{allPassGain, 441};
    allpass[3] = new AllPassDelay{allPassGain, 341};
    allpass[0] = new AllPassDelay{allPassGain, 225};
}

void MonoReverb::set(float room, float damp) {
    for (int i = 0; i < 8; i++) {
        combs[i]->room = room;
        if (damp < 1.0f) combs[i]->damp = damp;
    }
}

void MonoReverb::panic() {
    for (int i = 0; i < 8; i++) memset(combs[i]->buffer, 0, sizeof(float)*combSize);
    for (int i = 0; i < 8; i++) combs[i]->onepole = 0.0f;
    for (int i = 0; i < 4; i++) memset(allpass[i]->buffer, 0, sizeof(float)*frequency);
}

void MonoReverb::process(float & sample) {
    float out = 0;
    for (int j = 0; j < 8; j++) {
        out += combs[j]->process(sample);
    }
    sample = out;
    for (int i = 0; i < 4; i++) {
        allpass[i]->process(sample);
    }
}

StereoReverb::StereoReverb() {
    dry = 0.0f;
    for (int i = 0; i < 8; i++)
        r.combs[i]->n += 23;
    for (int i = 0; i < 4; i++) {
        r.allpass[i]->delay += 23;
    }
}

void StereoReverb::set(float room, float damp) {
    l.set(room, damp);
    r.set(room, damp);
}

void StereoReverb::process(float sample, StereoSample & stereoSample) {
    stereoSample.l = stereoSample.r = sample;
    r.process(stereoSample.r);
    stereoSample.r += dry * sample;
    l.process(stereoSample.l);
    stereoSample.l += dry * sample;
}

void Delay::panic() {
    std::fill(buffer, buffer + delayBufferSize, 0.0f);
    write = 0;
}

void Delay::process(float & sample) {
    float read = (float)write - (float)frequency * delay;
    if (read < 0) {
        read = (float)delayBufferSize + read;
    }
    float t = (float)fmod(read, 1);
    int read0 = (int)read;
    int read1 = ((int)read+1) % delayBufferSize;
    sample += buffer[read0] * t + buffer[read1] * (1 - t);
    buffer[write] = sample * decay;
    write++;
    write %= delayBufferSize;
}

void Fader::process(float sample, StereoSample & stereoSample) {
    float left = 0.5f + (balance * 0.5f);
    float right = 0.5f - (balance * 0.5f);
    stereoSample.l = left * sample;
    stereoSample.r = right * sample;
}

Sampler::Sampler(unsigned sampleCount) : sampleCount(sampleCount), overdubDecay(0.5f) {
    if (sampleCount > MaxSamples) throw std::runtime_error("sampler max exceeded");
    samples = new float[sampleCount]();
}

Sampler::~Sampler() {
    delete[] samples;
}

SamplerWriteHead::SamplerWriteHead(Sampler & sampler) : sampler(sampler), overdubDecay(1.0f), writeHead(0) { }

void SamplerWriteHead::write(float * samples, unsigned sampleCount) {
    unsigned readHead = 0;
    while (sampleCount) {
        sampler.samples[writeHead] *= overdubDecay; // naive attempt to prevent maxing out the sample and clipping
        sampler.samples[writeHead] += samples[readHead];
        readHead++;
        sampleCount--;
        writeHead++;
        if (writeHead >= sampler.sampleCount) writeHead = 0;
    }
}

SamplerReadHead::SamplerReadHead(Sampler & sampler) :speed(1.0f), sampler(sampler), readHead(0) {}

void SamplerReadHead::process(float * samples, unsigned sampleCount) {
    unsigned writeHead = 0;
    while (sampleCount) {
        samples[writeHead] = sampler.samples[readHead];
        writeHead++;
        readHead++;
        if (readHead > sampler.sampleCount) readHead = 0;
        sampleCount--;
    }
}
