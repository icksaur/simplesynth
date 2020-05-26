#include <iostream>
#include <cmath>
#include <ctime>
#include <SDL2/SDL.h>
#include <SDL2/SDL_main.h>

/* linker options: -lmingw32 -lSDLmain -lSDL -mwindows */

using namespace std;

const int noteA0 = 0;
const int noteA4 = 48;
//float notes[] = {27.50, 29.14, 30.87, 32.70, 34.65, 36.71, 38.89, 41.20, 43.65, 46.25, 49.00, 51.91, 55.00, 58.27, 61.74, 65.41, 69.30, 73.42, 77.78, 82.41, 87.31, 92.50, 98.00, 103.83, 110.00, 116.54, 123.47, 130.81, 138.59, 146.83, 155.56, 164.81, 174.61, 185.00, 196.00, 207.65, 220.00, 233.08, 246.94, 261.63, 277.18, 293.66, 311.13, 329.63, 349.23, 369.99, 392.00, 415.30, 440.00, 466.16, 493.88, 523.25, 554.37, 587.33, 622.25, 659.26, 698.46, 739.99, 783.99, 830.61, 880.00, 932.33, 987.77, 1046.50, 1108.73, 1174.66, 1244.51, 1318.51, 1396.91, 1479.98, 1567.98, 1661.22, 1760.00, 1864.66, 1975.53, 2093.00, 2217.46, 2349.32, 2489.02, 2637.02, 2793.83, 2959.96, 3135.96, 3322.44, 3520.00, 3729.31, 3951.07, 4186.01, 4434.92, 4698.64, 4978.03, 5274.04, 5587.65, 5919.91, 6271.93, 6644.88, 7040.00, };

const int frequency = 44100;
unsigned int sampleFrequency = 0;
unsigned int audioBufferSize = 0;
unsigned int outputAudioBufferSize = 0;

const unsigned short sampleRange = 32767;

int gate;

const float envelopePopThreshold = 0.003;


float noteFreqency(float a) {
    return 27.5f * pow(2.0f, (float)a);
}
float noteFreqency(int a, int halfSteps) {
    return noteFreqency((float)a + (float)halfSteps / 12.0f);
}

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
    void process(float * samples, int count, int gate) {
        float amp1;
        for (int i=0;i<count;i++) {
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
            samples[i] *= amp1;
            amp = amp1;
        }
    }
};

const float fade = 0.05;

struct Note {
    float glide;
    float note;
    float note0;
    float t;
    void generate(float * samples, int count) {
        if (note != note0 && t >= 1) {
            t = 0;
        }
        for (int i = 0; i < count; i++) {
            float n = (note * t + note0 * (1 - t));
            t += 1 / (float)frequency / glide;
            if (t >= 1) {
                t = 1;
                note0 = note;
            }
            samples[i] = n;
        }
    }
};

#define wrap(t) t = t > 1.0f ? 0.0f : t

struct Frequency {
    float t;
    void process(float * samples, int count) {
        for (int i = 0; i < count; i++) {
            samples[i] = noteFreqency(samples[i]);
            t += 1 / (float)samples[i];
        }
        wrap(t);
    }
};

SDL_Renderer * renderer;

int graph[400];

void renderState(float * samples, int count, float min = -1.0f, float max = 1.0f) {
    int width = count < 400 ? count : 400;
    float range = max - min;
    for (int i = 0; i < width; i++) {
        graph[i] = 400 - (samples[i] - min) / range * 400;
    }
}

// Generates waveforms which will bend each sample by sample frequency.
struct ModulatingSignal {
    float sine;
    float saw;
    float square;
    float noise;
    float t;
    void process(float * samples, int count) {
        for (int i=0; i<count; i++) {
            float h = sin(2.0 * M_PI * t) * sine;
            h += (2.0 * t - 1.0) * saw;
            h += (t < 0.5 ? -1 : 1) * square;
            h += ((float)rand() / (float)RAND_MAX * 2.0 - 1.0) * noise;
            t += samples[i] / frequency;
            samples[i] = h;
            t = fmod(t, 1.0f);
        }
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
    void process(float * samples, int count) {
        for (int i = 0; i < count; i++) {
            if (cursor >= wavelength - 1) {
                wavelength = frequency / samples[i];
                cursor = 0;
            }
            float t = (float)cursor / (float)(wavelength - 1);
            float h = sin(2.0 * M_PI * t) * sine;
            h += (2.0 * t - 1.0) * tri;
            h += (t < 0.5 ? -1 : 1) * square;
            h += ((float)rand() / (float)RAND_MAX * 2.0 - 1.0) * noise;
            samples[i] = h;
            cursor++;
        }
    }
};

struct LFO {
    float amplitude;
    float period;
    float t;
    void process(float * samples, int count) {
        for (int i = 0; i < count; i++) {
            samples[i] += sin(2.0 * M_PI * t) * amplitude;
            t += 1 / (float)frequency / period;
            wrap(t);
        }
    }
};

const int window = 7;

struct BasicFilter {
    float last[window];
    int next;
    float avg;
    float sum;

    void process(float * stream, int count) {
        for (int i = 0; i < count; i++) {
            sum += stream[i];
            sum -= last[next];
            last[next] = stream[i];
            stream[i] = sum / (float)window;
            ++next %= window;
        }
    }
};

const int delayBufferSize = 88211; // prime larger than 2 * 44100

#define nwrap(i) i < 0 ? (float)delayBufferSize - i : i

struct Delay {
    float delay;
    float decay;
    int write;
    float buffer[delayBufferSize];

    void process(float * samples, int count) {
        for (int i = 0; i < count; i++) {
            float read = (float)write - (float)frequency * delay;
            if (read < 0) {
                read = (float)delayBufferSize + read;
            }
            float t = fmod(read, 1);
            int read0 = (int)read;
            int read1 = ((int)read+1) % delayBufferSize;
            samples[i] += buffer[read0] * t + buffer[read1] * (1 - t);
            buffer[write] = samples[i] * decay;
            write++;
            write %= delayBufferSize;
        }
    }
};

// explicit delay in sample count
struct ExplicitDelay {
    int delay;
    float decay;
    int read;
    float buffer[delayBufferSize];
    void process(float * samples, int count) {
        delay = delay > (frequency - 1) ? (frequency - 1) : delay;
        for (int i = 0; i < count; i++) {
            samples[i] += buffer[read];
            buffer[(read + delay) % delayBufferSize] = samples[i] * decay;
            read++;
            read %= delayBufferSize;
        }
    }
};

int baseDelay = 441;

struct Reverb {
    ExplicitDelay delays[10];
    Reverb() {
        delays[0].delay = baseDelay + 29; delays[0].decay = 0.1;
        delays[1].delay = baseDelay + 601; delays[1].decay = 0.2;
        delays[2].delay = baseDelay + 1291; delays[2].decay = 0.4;
        delays[3].delay = baseDelay + 2053; delays[3].decay = 0.6;
        delays[4].delay = baseDelay + 2819; delays[4].decay = 0.7;
        delays[5].delay = baseDelay + 3643; delays[5].decay = 0.8;
        delays[6].delay = baseDelay + 4493; delays[6].decay = 0.7;
        delays[7].delay = baseDelay + 5387; delays[7].decay = 0.6;
        delays[8].delay = baseDelay + 6221; delays[8].decay = 0.2;
        delays[9].delay = baseDelay + 7103; delays[9].decay = 0.1;
    }

    void process(float * samples, int count) {
        for (int d = 0; d < 4; d++) {
            delays[d].process(samples, count);
        }
    }
};

struct CombinedReverb {
    float decay;
    int delays[10];
    int read;
    float buffer[delayBufferSize];
    CombinedReverb() {
        decay = 0.1;
        delays[0] = 29;
        delays[0] = 601;
        delays[1] = 1291;
        delays[2] = 2053;
        delays[3] = 2819;
        delays[4] = 3643;
        delays[5] = 4493;
        delays[6] = 5387;
        delays[7] = 6221;
        delays[8] = 7103;
    }

    void process(float * samples, int count) {
       for (int i = 0; i < count; i++) {
            samples[i] += buffer[read];
            for (int j = 4; j < 9; j++) {
                buffer[(read + delays[j]) % delayBufferSize] += samples[i] * decay;
            }
            read++;
            read %= delayBufferSize;
        }
    }
};

Note noteGen{ 0.001, 4, 4 };
Frequency noteFreq;
ModulatingSignal modsig { 0, 1, 0, 0 };
FixedSignal fixedsig { 0, 1, 0, 0 };
Envelope env{0, 0, 1, 0, 0, 0, 4};
BasicFilter basicFilter;
Delay delay{ 0.5, 0.5 };
Reverb reverb;
CombinedReverb combinedReverb;

ExplicitDelay ed1{ 229, 0.5 };
ExplicitDelay ed2{ 1069, 0.5 };
ExplicitDelay ed3{ 3181, 0.5 };
ExplicitDelay ed4{ 6053, 0.5 };
ExplicitDelay ed5{ 7919, 0.5 };
LFO lfo{ 0.01, 0.5 };

int mode = 0;

#define FLOATSTREAM (float*)stream, len / 4

void fillAudio(void *unused, Uint8 *stream, int len) {
    noteGen.generate(FLOATSTREAM);
    lfo.process(FLOATSTREAM);
    //renderState((float*)stream, len / 4, noteGen.note-0.1, noteGen.note+0.1);
    noteFreq.process(FLOATSTREAM);
    modsig.process(FLOATSTREAM);
    env.process(FLOATSTREAM, gate);
    if (mode == 0) {
        combinedReverb.process(FLOATSTREAM);
    }
    else if (mode == 1) {
        reverb.process(FLOATSTREAM);
    }
    basicFilter.process(FLOATSTREAM);
    renderState((float*)stream, len);
}

void noteOn(SDL_Event & e, int note) {
    if (e.key.repeat) return;
    gate++;
    //freqGen.note = (float)note / 12.0f;
    noteGen.note = (float)note / 12.0f;

    env.trigger();
}

int main(int argc, char *argv[]) {
    srand(time(0));
    if( SDL_Init(SDL_INIT_TIMER | SDL_INIT_AUDIO ) <0 ) {
        cout << "Unable to init SDL: " << SDL_GetError() << endl;
        return 1;
    }

    int format = AUDIO_F32;

    SDL_AudioSpec desired, obtained;
    desired.freq = 44100;
    desired.format = format;
    desired.samples = 1024;
    desired.callback = fillAudio;
    desired.userdata = NULL;
    desired.channels = 1;

    if ( SDL_OpenAudio(&desired, &obtained) < 0 ) {
        fprintf(stderr, "AudioMixer, Unable to open audio: %s\n", SDL_GetError());
        return 1;
    }

    audioBufferSize = obtained.samples;
    sampleFrequency = obtained.freq;

    /* if the format is 16 bit, two bytes are written for every sample */
    if (obtained.format != format) {
        fprintf(stderr, "did not get desired format\n");
        return 1;
    }

    SDL_Window * window = SDL_CreateWindow(
        "audio", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, 400, 400, 0);
    renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);

    SDL_PauseAudio(0);

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
                case SDLK_w: noteOn(event, octave+1); break;
                case SDLK_s: noteOn(event, octave+2); break;
                case SDLK_e: noteOn(event, octave+3); break;
                case SDLK_d: noteOn(event, octave+4); break;
                case SDLK_f: noteOn(event, octave+5); break;
                case SDLK_t: noteOn(event, octave+6); break;
                case SDLK_g: noteOn(event, octave+7); break;
                case SDLK_y: noteOn(event, octave+8); break;
                case SDLK_h: noteOn(event, octave+9); break;
                case SDLK_u: noteOn(event, octave+10); break;
                case SDLK_j: noteOn(event, octave+11); break;
                case SDLK_k: noteOn(event, octave+12); break;
                case SDLK_F4: if (event.key.keysym.mod & KMOD_ALT) running = false; break;
                case SDLK_PAGEDOWN: octave = octave ? octave-12 : octave; break;
                case SDLK_PAGEUP: octave = octave < (12*7) ? octave+12 : octave; break;
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
            }
        }
        SDL_Delay(10);
        SDL_SetRenderDrawColor(renderer, 0, 0, 0, 0);
        SDL_RenderClear(renderer);
        SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
        for (int i = 0; i < 400; i++)
            SDL_RenderDrawLine(renderer, i, 400, i, graph[i]);
        SDL_RenderPresent(renderer);
    }
    printf("quitting\n");
    SDL_Quit();
    return 0;
}
