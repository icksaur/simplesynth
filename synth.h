#pragma once

#include <cmath>
#include <stdexcept>
#include <memory>
#include <cstring>

const int noteA4 = 48;
const int frequency = 44100;
const float envelopePopThreshold = 0.003f;
const int delayBufferSize = 88211; // prime larger than 2 * 44100
const float defaultRoom = 0.84f;
const float allPassGain = 0.618f; // reciporical of golden ratio
const unsigned MaxSamples = frequency * 60;

float noteFrequency(float a);
float noteFrequency(int a, int halfSteps);

// Represents a stereo audio sample with left and right channels
struct StereoSample { float l, r; };

// ADSR envelope generator for controlling amplitude over time
struct Envelope {
    float a, d, s, r;
    float t;
    int _gate;
    char state; 
    float amp;
    void trigger();
    float generate(int gate);
};

// Note generator with glide/portamento between pitches
struct Note {
    float glide;
    float note;
    float note0;
    float t;
    float generate();
};

#define wrap(t) t = t > 1.0f ? t-1.0f : t

// Generates waveforms which will bend each sample by sample frequency.
struct Oscillator {
    float sine;
    float saw;
    float square;
    float noise;
    float t;
    float generate(float f);
};

// Low frequency oscillator for modulation effects
struct LFO {
    float amplitude;
    float period;
    float t;
    float generate();
};

#define MIN(a, b) (a < b ? a : b)

// The Scientist and Engineer's Guide to Digital Signal Processing, ch19, Steven W Smith
// cutoff frequency should be f=0.5..1.0 : e^-2PIf (exp(2.0f * M_PI * f))
struct StagedFilter {
    float k; // if this goes to 1.0, y becomes stuck
    float y[4];
    StagedFilter();
    void process(float & sample);
};

typedef StagedFilter LowpassFilter;

#undef MIN

// Bandpass filter for isolating specific frequency ranges
struct BandpassFilter {
    float f; // frequency (resonant frequency / 44100)
    float r; // bandwidth
    float x[2]; // input history
    float y[2]; // output history
    BandpassFilter();
    void process(float & sample);
};

// Classic Moog-style low/high/bandpass filter with resonance
struct MoogFilter {
    float frequency, resonance;
    enum FilterType { low, high, band } type;
    float f, p, q;
    float b0, b1, b2, b3, b4;
    float t1, t2;
    MoogFilter(float frequency, float resonance);
    void panic();
    void set(float frequency, float resonance);
    void process(float & sample);
};

float MoogFreq(float f);

typedef BandpassFilter FormantBand;

float peak(int hz);
float bandwidth(int hz);
float gain(int db);

// Defines formant frequencies for vowel sounds with frequency, amplitude, and bandwidth
struct FormantVowel {
    int hz[5];
    int db[5];
    int bw[5];
};

// External declarations for formant vowels
extern FormantVowel soprano_a, soprano_e, soprano_i, soprano_o, soprano_u;
extern FormantVowel alto_a, alto_e, alto_i, alto_o, alto_u;
extern FormantVowel countertenor_a, countertenor_e, countertenor_i, countertenor_o, countertenor_u;
extern FormantVowel tenor_a, tenor_e, tenor_i, tenor_o, tenor_u;
extern FormantVowel bass_a, bass_e, bass_i, bass_o, bass_u;

// Vocal filter that simulates human vowel sounds using formant filtering
struct VocalFilter {
    int formants;
    FormantBand filters[5];
    VocalFilter();
    void setVowelFormants(FormantVowel * vowel);
    void process(float & sample);
};

// All-pass filter with delay for creating reverb effects
struct AllPassDelay {
    float k;
    int delay;
    float buffer[frequency];
    int read;
    void process(float & sample);
};

// A building block for other things.
struct AllPassFilter {
    float k;
    float r;
    void process(float & sample);
};

// from explanation at https://www.dsprelated.com/freebooks/pasp/Freeverb.html
const int combSize = 3000;
struct LowpassFeedbackCombFilter {
    float damp, room;
    int n;
    int read;
    float * buffer;
    float onepole;
    LowpassFeedbackCombFilter(float room, float damp, int delay);
    ~LowpassFeedbackCombFilter();
    float process(float x);
};

typedef LowpassFeedbackCombFilter LBCF;

// Mono reverb engine using comb and all-pass filters
struct MonoReverb {
    LBCF * combs[8];
    AllPassDelay * allpass[4];
    MonoReverb();
    void set(float room, float damp);
    void panic();
    void process(float & sample);
};

// Stereo reverb with independent left/right processing and dry mix
struct StereoReverb {
    MonoReverb l;
    MonoReverb r;
    float dry;
    StereoReverb();
    void set(float room, float damp);
    void process(float sample, StereoSample & stereoSample);
};

#define nwrap(i) i < 0 ? (float)delayBufferSize - i : i

// Echo/delay effect with configurable delay time and feedback
struct Delay {
    float delay;
    float decay;
    int write;
    float buffer[delayBufferSize];

    void panic();
    void process(float & sample);
};

// Stereo balance/panning control for positioning audio in stereo field
struct Fader {
    float balance;
    void process(float sample, StereoSample & stereoSample);
};

// Audio sample buffer for recording and playback functionality
struct Sampler {
    float * samples;
    float overdubDecay;
    unsigned sampleCount;
    Sampler() = delete;
    Sampler(unsigned sampleCount);
    ~Sampler();
};

// Write head for recording audio into a sampler buffer
struct SamplerWriteHead {
    Sampler & sampler;
    float overdubDecay;
    unsigned writeHead;
    SamplerWriteHead(Sampler & sampler);
    void write(float * samples, unsigned sampleCount);
};

// Read head for playing back audio from a sampler buffer
struct SamplerReadHead {
    Sampler & sampler;
    float speed;
    unsigned readHead;
    SamplerReadHead(Sampler & sampler);
    void process(float * samples, unsigned sampleCount);
};