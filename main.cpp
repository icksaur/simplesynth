#include <cmath>
#include <ctime>
#include <SDL2/SDL.h>
#include <SDL2/SDL_main.h>

#include "GLObject.h"
#include "glxw.h"
#include "gui.h"

using namespace std;

const int noteA4 = 48;
const int frequency = 44100;
unsigned int sampleFrequency = 0;
unsigned int audioBufferSize = 0;

int gate;

const float envelopePopThreshold = 0.003;
const int delayBufferSize = 88211; // prime larger than 2 * 44100

const int rez = 512;

SDL_Renderer * renderer;
int graph[rez];
bool graphDirty = false;
void renderState(float * samples, int count) {
    const float min = -1.0f;
    const float max = 1.0f;
    const float range = max - min;
    if (!graphDirty) return;
    graphDirty = false;
    int width = count < rez ? count : rez;
    for (int i = 0; i < width; i++) {
        graph[i] = rez - (samples[i] - min) / range * rez;
    }
}

float noteFreqency(float a) {
    return 27.5f * pow(2.0f, (float)a);
}
float noteFreqency(int a, int halfSteps) {
    return noteFreqency((float)a + (float)halfSteps / 12.0f);
}

struct StereoSample { float l, r; };

struct Envelope {
    float a, d, s, r;
    float t;
    int _gate;
    char state; 
    float amp;
    void trigger() {
        _gate = 1;
        state = 0; 
        t = 0;
    }
    void process(float & sample, int gate) {
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
        if (abs(ampDif) > envelopePopThreshold) {
            amp1 = amp + copysign(envelopePopThreshold, ampDif);
        }
        sample *= amp1;
        amp = amp1;
    }
};

struct Note {
    float glide;
    float note;
    float note0;
    float t;
    void generate(float * samples, int count) {
        for (int i = 0; i < count; i++) {
            generate(samples[i]);
        }
    }
    void generate(float & sample) {
        float n = (note * t + note0 * (1 - t));
        t += 1 / (float)frequency / glide;
        if (t >= 1) {
            t = 1;
            note0 = note;
        }
        sample = n;
    }
};

//#define wrap(t) t = t > 1.0f ? 0.0f : t
#define wrap(t) t = t > 1.0f ? t-1.0f : t

struct Frequency {
    float t;
    void process(float & sample) {
        sample = noteFreqency(sample);
        t += 1 / sample;
        wrap(t);
    }
};

// Generates waveforms which will bend each sample by sample frequency.
struct ModulatingSignal {
    float sine;
    float saw;
    float square;
    float noise;
    float t;
    void process(float & sample) {
        float h = sin(2.0f * M_PI * t) * sine;
        h = t / 0.5f - 1.0f;
        h += (t < 0.5 ? -1 : 1) * square;
        h += ((float)rand() / (float)RAND_MAX * 2.0 - 1.0) * noise;
        t += sample / frequency;
        sample = h;
        wrap(t);
    }
};


// Waveform will not change until period is over. 
// Period is snapped to sample count.
struct FixedSignal {
    float sine;
    float tri;
    float square;
    float noise;
    int wavelength;
    int cursor;
    void process(float & sample) {
        if (cursor >= wavelength - 1) {
            wavelength = frequency / sample;
            cursor = 0;
        }
        float t = (float)cursor / (float)(wavelength - 1);
        float h = sin(2.0 * M_PI * t) * sine;
        h += (2.0 * t - 1.0) * tri;
        h += (t < 0.5 ? -1 : 1) * square;
        h += ((float)rand() / (float)RAND_MAX * 2.0 - 1.0) * noise;
        sample = h;
        cursor++;
    }
};

struct LFO {
    float amplitude;
    float period;
    float t;
    void process(float & sample) {
        sample += sin(2.0 * M_PI * t) * amplitude;
        t += 1 / (float)frequency / period;
        wrap(t);
    }
};

const int window = 7;

// A delay with no feedback (first two statements of process are reversed)
struct CombFilter {
    float decay;
    int delay;
    int read;
    float buffer[delayBufferSize];
    float process(float sample) {
        buffer[(read + delay) % delayBufferSize] = sample * decay;
        sample += buffer[read];
        read++;
        read %= delayBufferSize;
        return sample;
    }
};

struct BasicFilter {
    float k;
    float r;
    void process(float & sample) {
        sample = (sample - r) * k + r;
        r = sample;
    }
};

// The Scientist and Engineer's Guide to Digital Signal Processing, ch19, Steven W Smith
// cutoff frequency should be f=0.5..1.0 : e^-2PIf (exp(2.0f * M_PI * f))
struct StagedFilter {
    float k; // if this goes to 1.0, y becomes stuck
    float y[4];
    StagedFilter() {
        k = 0.5f;
    }
    void process(float & sample) {
        float k = min(this->k, 0.95f); // cache local K less than 0.5f
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
};

// The Scientist and Engineer's Guide to Digital Signal Processing, ch19
struct BandpassFilter {
    float k; // (0.0 to 0.5] if it goes to zero, it gets stuck
    float k0; // previous k
    float BW; // bandwidth
    float x[2]; // input history
    float y[2]; // output history
    float R;
    float K;
    float co; // unnamed coefficient which is common
    void coeff() { // recalculate K, R, and co
        if (k == k0) return;
        k = max(k, 0.0001f); // prevent filter getting stuck
        R = 1.0f - 3.0f * BW;
        co = 2.0f * cos(2.0f * M_PI * k);
        K = (1.0f - R * co + (R * R)) / (2.0f - co);
        k0 = k;
    }
    BandpassFilter() {
        k = 0.2f;
        BW = 0.005f;
    }
    void process(float & sample) {
        coeff();
        float x0 = sample;
        sample = sample * (1.0f - K);
        sample += x[0] * (K - R) * co;
        sample += x[1] * (R * R - K);
        sample += y[0] * R * co;
        sample += y[1] * -R * R;
        
        y[1] = y[0];
        y[0] = sample;
        x[1] = x[0];
        x[0] = x0;
    }
};

struct BandpassFilter2 {
    float f; // frequency (resonant frequency / 44100)
    float r; // bandwidth
    float x[2]; // input history
    float y[2]; // output history
    BandpassFilter2() {
        f = 0.2f;
        r = 0.95f;
    }
    void process(float & sample) {
        float x0 = sample;
        float b0 = (1.0f - r * r) / 2.0f;
        sample = sample * b0;
        sample += x[1] * -b0;
        sample -= y[0] * -2.0f * r * cos(2.0f * M_PI * f);
        sample -= y[1] * r * r;
        
        y[1] = y[0];
        y[0] = sample;
        x[1] = x[0];
        x[0] = x0;
    }
};

struct MoogFilter {
    float frequency, resonance;
    enum FilterType { low, high, band } type;
    float f, p, q;
    float b0, b1, b2, b3, b4;
    float t1, t2;
    MoogFilter(float _frequency, float _resonance) : frequency(_frequency),resonance(_resonance) {
        set();
    }
    void set() {
        q = 1.0f - frequency;
        p = frequency + 0.8f * frequency * q;
        f = p + p - 1.0f;
        q = resonance * (1.0f + 0.5f * q * (1.0f - q + 5.6f * q * q));
    }
    void process(float & sample) {
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
};

float MoogFreq(float f) {
    return exp(2.0f * M_PI * f);
}

typedef BandpassFilter2 FormantBand;

float peak(int hz) {
    return (float)hz / (float)frequency;
}
float bandwidth(int hz) {
    return 1.0f - (float)hz / (float)frequency;
}

struct VocalFilter {
    float dry;
    FormantBand filters[3];
    //VocalFilter()  : filters{ FormantBand(0.1, 0.9), FormantBand(0.2, 0.9), FormantBand(0.3, 0.9) } {
    VocalFilter() {
        /*
        filters[0].k = 730.0f / (float)frequency;
        filters[1].k = 1090.0f / (float)frequency;
        filters[2].k = 2240.0f / (float)frequency;
        */
        filters[0].f = peak(800); filters[0].r = bandwidth(80);
        filters[1].f = peak(1150); filters[1].r = bandwidth(90);
        filters[2].f = peak(2900); filters[2].r = bandwidth(120);
        /*
        filters[0].set();
        filters[1].set();
        filters[2].set();
        */
    }
    void process(float & sample) {
        float x0 = sample, x = sample;
        sample *= dry;
        filters[0].process(x); sample += x; x = x0;
        filters[1].process(x); sample += x; x = x0;
        filters[2].process(x); sample += x; x = x0;
    }
};

typedef StagedFilter LowpassFilter;

// explicit delay in sample count
struct ExplicitDelay {
    int delay;
    float decay;
    int read;
    float buffer[delayBufferSize];
    void process(float & sample) {
        sample += buffer[read];
        buffer[(read + delay) % delayBufferSize] = sample * decay;
        read++;
        read %= delayBufferSize;
    }
};

struct AllPassDelay {
    float k;
    int delay;
    float buffer[frequency];
    int read;
    void process(float & sample) {
        float r1 = buffer[read] + sample * k; // peek in delay
        sample = r1 * -k + buffer[read];
        buffer[(read + delay) % frequency] = r1 * 0.5;
        ++read %= frequency;
    }
};

// Not sure what this does. It's a building block for other things.
struct AllPassFilter {
    float k;
    float r;
    void process(float & sample) {
        float r1 = sample + r * k;
        sample = r1 * -k + r;
        r = r1;
    }
};

// from explanation at https://www.dsprelated.com/freebooks/pasp/Freeverb.html
const int combSize = 3000;
struct LowpassFeedbackCombFilter {
    float damp, room;
    int n;
    int read;
    float * buffer;
    float onepole;
    LowpassFeedbackCombFilter(float room, float damp, int delay) {
        if (delay >= combSize) throw std::runtime_error("delay is too long");
        this->damp = damp;
        this->room = room;
        this->n = delay;
        buffer = new float[combSize]; // guess this could be n
        memset(buffer, 0, sizeof(float) * combSize);
        read = combSize - n;
        onepole = 0.0f;
    }
    ~LowpassFeedbackCombFilter() {
        delete[] buffer;
    }
    float process(float x) {
        float y = buffer[read]; // delay line
        onepole = (1.0f - damp) / (1.0f - damp * onepole);
        int write = (read + n) % combSize;
        buffer[write] = x + room * onepole * y;
        read = (read + 1) % combSize;
        return y;
    }
};

typedef LowpassFeedbackCombFilter LBCF;

const float defaultRoom = 0.84f;
const float allPassGain = 0.618f; // reciporical of golden ratio

struct MonoReverb {
    LBCF * combs[8];
    AllPassDelay * allpass[4];
    MonoReverb() {
        combs[0] = new LBCF(defaultRoom, .2, 1617);
        combs[1] = new LBCF(defaultRoom, .2, 1557);
        combs[2] = new LBCF(defaultRoom, .2, 1491);
        combs[3] = new LBCF(defaultRoom, .2, 1277);
        combs[4] = new LBCF(defaultRoom, .2, 1422);
        combs[5] = new LBCF(defaultRoom, .2, 1356);
        combs[6] = new LBCF(defaultRoom, .2, 1188);
        combs[7] = new LBCF(defaultRoom, .2, 1116);
        allpass[1] = new AllPassDelay{allPassGain, 556};
        allpass[2] = new AllPassDelay{allPassGain, 441};
        allpass[3] = new AllPassDelay{allPassGain, 341};
        allpass[0] = new AllPassDelay{allPassGain, 225};
    }
    void set(float room, float damp) {
        for (int i = 0; i < 8; i++) {
            combs[i]->room = room;
            if (damp < 1.0f) combs[i]->damp = damp;
        }
    }
    void panic() {
        for (int i = 0; i < 8; i++) memset(combs[i]->buffer, 0, sizeof(float)*combSize);
        for (int i = 0; i < 8; i++) combs[i]->onepole = 0.0f;
        for (int i = 0; i < 4; i++) memset(allpass[i]->buffer, 0, sizeof(float)*frequency);
    }
    void process(float & sample) {
        float out = 0;
        for (int j = 0; j < 8; j++) {
            out += combs[j]->process(sample);
        }
        sample = out;
        for (int i = 0; i < 4; i++) {
            allpass[i]->process(sample);
        }
    }
};

struct StereoReverb {
    MonoReverb l;
    MonoReverb r;
    float dry;
    StereoReverb() {
        dry = 0.0f;
        for (int i = 0; i < 8; i++)
            r.combs[i]->n += 23;
        for (int i = 0; i < 4; i++) {
            r.allpass[i]->delay += 23;
        }
    }
    void process(float sample, StereoSample & stereoSample) {
        stereoSample.l = stereoSample.r = sample;
        r.process(stereoSample.r);
        stereoSample.r += dry * sample;
        l.process(stereoSample.l);
        stereoSample.l += dry * sample;
    }
};

#define nwrap(i) i < 0 ? (float)delayBufferSize - i : i

struct Delay {
    float delay;
    float decay;
    int write;
    float buffer[delayBufferSize];

    void process(float & sample) {
        float read = (float)write - (float)frequency * delay;
        if (read < 0) {
            read = (float)delayBufferSize + read;
        }
        float t = fmod(read, 1);
        int read0 = (int)read;
        int read1 = ((int)read+1) % delayBufferSize;
        sample += buffer[read0] * t + buffer[read1] * (1 - t);
        buffer[write] = sample * decay;
        write++;
        write %= delayBufferSize;
    }
};

struct Fader {
    float balance;
    void process(float sample, StereoSample & stereoSample) {
        float left = 0.5f + (balance * 0.5f);
        float right = 0.5f - (balance * 0.5f);
        stereoSample.l = left * sample;
        stereoSample.r = right * sample;
    }
};

int baseDelay = 441;

Note noteGen{ 0.001, 4, 4 };
Frequency noteFreq;
ModulatingSignal modsig { 0, 1, 0, 0 };
FixedSignal fixedsig { 0, 1, 0, 0 };
Envelope env{0.05f, 0, 1, 0.05f, 0, 0, 4};
LowpassFilter lowpassFilter;
AllPassFilter allPassFilter{ 0.5 };
AllPassDelay allPassDelay{ -0.5, 22050 };
Delay delay{ 0.5, 0.5 };
ExplicitDelay ed1{ 229, 0.5 };
ExplicitDelay ed2{ 1069, 0.5 };
ExplicitDelay ed3{ 3181, 0.5 };
ExplicitDelay ed4{ 6053, 0.5 };
ExplicitDelay ed5{ 7919, 0.5 };
StereoReverb stereoReverb;
BandpassFilter bandpass;
MoogFilter moogFilter(0.5f, 0.0f);
VocalFilter vocalFilter;
LFO lfo{ 0.01, 0.5 };
Fader fader;

int mode = 0;

float monoSamples[1024];
#define FLOATSTREAM (float*)monoSamples, len / 4

void fillAudio(void *unused, Uint8 *stream, int len) {
    float * floatStream = (float*)stream;
    StereoSample * stereoSamples = (StereoSample*)stream;
    for (unsigned i = 0; i < len/sizeof(*stereoSamples); i++) {
        float sample;
        noteGen.generate(sample);
        lfo.process(sample);
        noteFreq.process(sample);
        modsig.process(sample);
        env.process(sample, gate);
        if (mode == 0) {
            vocalFilter.process(sample);
            lowpassFilter.process(sample);
            fader.process(sample, stereoSamples[i]);
        } else if (mode == 1) {
            vocalFilter.process(sample);
            stereoReverb.process(sample, stereoSamples[i]);
            lowpassFilter.process(stereoSamples[i].r);
            lowpassFilter.process(stereoSamples[i].l);
        } else {
            //moogFilter.process(sample);
            fader.process(sample, stereoSamples[i]);
        }
    }
    renderState((float*)stream, len);
}

void noteOn(SDL_Event & e, int note) {
    if (e.key.repeat) return;
    gate++;
    noteGen.note = (float)note / 12.0f;
    env.trigger();
}

#define WHITE 255,255,255,255
#define BLACK 0,0,0,255
const char * fullScreenTexturedQuadSource = R"(
#if VERTEX_SHADER
layout(location=0) in vec2 vertex;
layout(location=1) in vec2 uv;
out vec2 uvco;
void main(void) {
uvco = uv;
gl_Position = vec4(vertex.xy, 0, 1);
}
#elif FRAGMENT_SHADER
layout(binding = 0) uniform sampler2DArray textures;
in vec2 uvco;
void main() {
    gl_FragColor = texture(textures, vec3(uvco.st, 0));
}
#endif
)";

int step = 20;
int elementY;
bool param(GUI & gui, float & value, const char * label) {
    bool result = gui.slider(value, 0, elementY, 100, 20);
    gui.label(110, elementY + 10, label);
    elementY += step;
    return result;
}

int main(int argc, char *argv[]) {
    srand(time(0));
    int xrez = rez, yrez = rez;
    SDL_Init(SDL_INIT_VIDEO | SDL_INIT_AUDIO);
    SDL_Window * window = SDL_CreateWindow("test", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, xrez, xrez, SDL_WINDOW_OPENGL | SDL_WINDOW_ALLOW_HIGHDPI);
    SDL_GLContext context = SDL_GL_CreateContext(window);

    struct OnExit {
        SDL_Window * window;
        SDL_GLContext context;
        ~OnExit() {
            SDL_GL_DeleteContext(context);
            SDL_DestroyWindow(window);
        }
    } cleanup{ window, context };

    if (glxwInit()) return 0;

    int format = AUDIO_F32;
    SDL_AudioSpec desired, obtained;
    desired.freq = 44100;
    desired.format = format;
    desired.samples = 1024;
    desired.callback = fillAudio;
    desired.userdata = NULL;
    desired.channels = 2;
    if ( SDL_OpenAudio(&desired, &obtained) < 0 ) {
        fprintf(stderr, "AudioMixer, Unable to open audio: %s\n", SDL_GetError());
        return 1;
    }
    audioBufferSize = obtained.samples;
    sampleFrequency = obtained.freq;
    if (obtained.format != format) {
        fprintf(stderr, "did not get desired format\n");
        return 1;
    }
    SDL_PauseAudio(0);

    ShaderProgram program(fullScreenTexturedQuadSource);
    VAOBuilder builder(1);
    builder.addFloats(0, 0, 2).addFloats(0, 1, 2);
    VAO vao(program, builder);
    const float w = 1.0f;
    float vertices[] = { -w,-w, 0,0, -w,w, 0,1, w,w, 1,1, w,-w, 1,0 };
    vao.vertexData(vertices, sizeof(vertices));

    unsigned graphBytesSize = rez * rez * 4;
    char * graphBytes = new char[graphBytesSize];
    unsigned * graphTexels = (unsigned*)graphBytes;
    memset(graphBytes, 0, graphBytesSize);
    ArrayTexture graphTexture(rez, 1, graphBytes, graphBytesSize);

    GUI gui;
    int octave = noteA4;
    bool running = true;

    SDL_Event event;
    while (running) {
        while (SDL_PollEvent(&event)) {
            switch (event.type) {
            case SDL_KEYDOWN:
                switch (event.key.keysym.sym) {
                case SDLK_ESCAPE: running = false; break;
                case SDLK_a: noteOn(event, octave); break;
                case SDLK_w: noteOn(event, octave + 1); break;
                case SDLK_s: noteOn(event, octave + 2); break;
                case SDLK_e: noteOn(event, octave + 3); break;
                case SDLK_d: noteOn(event, octave + 4); break;
                case SDLK_f: noteOn(event, octave + 5); break;
                case SDLK_t: noteOn(event, octave + 6); break;
                case SDLK_g: noteOn(event, octave + 7); break;
                case SDLK_y: noteOn(event, octave + 8); break;
                case SDLK_h: noteOn(event, octave + 9); break;
                case SDLK_u: noteOn(event, octave + 10); break;
                case SDLK_j: noteOn(event, octave + 11); break;
                case SDLK_k: noteOn(event, octave + 12); break;
                case SDLK_F4: if (event.key.keysym.mod & KMOD_ALT) running = false; break;
                case SDLK_PAGEDOWN: octave = octave ? octave - 12 : octave; break;
                case SDLK_PAGEUP: octave = octave < (12 * 7) ? octave + 12 : octave; break;
                case SDLK_SPACE: ++mode %= 3; break;
                }
                break;
            case SDL_KEYUP:
                switch (event.key.keysym.sym) {
                case SDLK_a: gate--; break;
                case SDLK_w: gate--; break;
                case SDLK_s: gate--; break;
                case SDLK_e: gate--; break;
                case SDLK_d: gate--; break;
                case SDLK_f: gate--; break;
                case SDLK_t: gate--; break;
                case SDLK_g: gate--; break;
                case SDLK_y: gate--; break;
                case SDLK_h: gate--; break;
                case SDLK_u: gate--; break;
                case SDLK_j: gate--; break;
                case SDLK_k: gate--; break;
                }
                break;
            case SDL_QUIT:
                running = false;
                break;
            case SDL_MOUSEMOTION: gui.pushMouseMotion(event.motion.x, rez - event.motion.y); break;
            case SDL_MOUSEBUTTONDOWN: if (event.button.button == SDL_BUTTON_LEFT) gui.pushMouseDown(0, event.button.x, rez - event.button.y); break;
            case SDL_MOUSEBUTTONUP: if (event.button.button == SDL_BUTTON_LEFT) gui.pushMouseUp(0, event.button.x, rez - event.button.y); break;
            }
        }

        memset(graphTexels, 0, sizeof(unsigned) * rez * rez);
        for (int i = 0; i < rez; i++) {
            for (int j = 0; j < rez; j++) {
                if (graph[i] < j) {
                    graphTexels[rez * j + i] = 0xFFFFFFFF;
                }
            }
        }
        glActiveTexture(GL_TEXTURE0 + 0);
        graphTexture.bind();
        glTexSubImage3D(GL_TEXTURE_2D_ARRAY, 0, 0, 0, 0, rez, rez, 1, GL_RGBA, GL_UNSIGNED_BYTE, graphTexels);

        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        program.use();
        vao.bind();
        glDrawArrays(GL_TRIANGLE_FAN, 0, 4);

        static int x, y; // window drag
        static char krtString[12];
        gui.start(x, y, 250, 300, "params");
        elementY = 0;
        static float room = 0.82f;
        static float damp = 0.2f;
        static float fade = 0.5f;

        param(gui, modsig.sine, "sine");
        param(gui, modsig.square, "square");
        param(gui, modsig.saw, "saw");
        param(gui, modsig.noise, "noise");

        static float f1, f2, f3 = 0.1f;

        const char * filterNames[3]{ "low", "high", "band" };

        switch (mode) {
        case 0:
            param(gui, vocalFilter.dry, "dry");
            /*
            if (param(gui, vocalFilter.filters[0].frequency, "formant1")) vocalFilter.filters[0].set();
            if (param(gui, vocalFilter.filters[1].frequency, "formant2")) vocalFilter.filters[1].set();
            if (param(gui, vocalFilter.filters[2].frequency, "formant3")) vocalFilter.filters[2].set();
            */
            if (param(gui, f1, "formant1")) vocalFilter.filters[0].f = f1;
            if (param(gui, f2, "formant2")) vocalFilter.filters[1].f = f2;
            if (param(gui, f3, "formant3")) vocalFilter.filters[2].f = f3;

            if (param(gui, fade, "fader")) fader.balance = (fade - 0.5f) * -2.0f;
            break;
        case 1:
            param(gui, lowpassFilter.k, "lowpass");
            if (param(gui, room, "room")) {
                stereoReverb.l.set(room, damp);
                stereoReverb.r.set(room, damp);
            }
            if (param(gui, damp, "damp")) {
                stereoReverb.l.set(room, damp);
                stereoReverb.r.set(room, damp);
            }
            param(gui, stereoReverb.dry, "dry");
            if (gui.button(0, elementY, 100, 20, "panic")) {
                stereoReverb.l.panic();
                stereoReverb.r.panic();
            }
            break;
        case 2:
            if (param(gui, moogFilter.frequency, "filter")) {
                moogFilter.set();
            }
            if (param(gui, moogFilter.resonance, "resonance")) {
                moogFilter.set();
            }
            if (gui.button(0, elementY, 100, 20, filterNames[moogFilter.type])) {
                moogFilter.type = (MoogFilter::FilterType)((moogFilter.type + 1) % 3);
            }
        }
        gui.end();
        gui.render();
        graphDirty = true;
        SDL_GL_SwapWindow(window);
        SDL_Delay(1000/60);
    }

    SDL_GL_DeleteContext(context);
    SDL_DestroyWindow(window);
    return 0;
}
