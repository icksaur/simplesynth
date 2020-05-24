#include <iostream>
#include <cmath>
#include <SDL2/SDL.h>
#include <SDL2/SDL_main.h>

/* linker options: -lmingw32 -lSDLmain -lSDL -mwindows */

using namespace std;

const int noteA0 = 0;
const int noteA4 = 48;
float notes[] = {27.50, 29.14, 30.87, 32.70, 34.65, 36.71, 38.89, 41.20, 43.65, 46.25, 49.00, 51.91, 55.00, 58.27, 61.74, 65.41, 69.30, 73.42, 77.78, 82.41, 87.31, 92.50, 98.00, 103.83, 110.00, 116.54, 123.47, 130.81, 138.59, 146.83, 155.56, 164.81, 174.61, 185.00, 196.00, 207.65, 220.00, 233.08, 246.94, 261.63, 277.18, 293.66, 311.13, 329.63, 349.23, 369.99, 392.00, 415.30, 440.00, 466.16, 493.88, 523.25, 554.37, 587.33, 622.25, 659.26, 698.46, 739.99, 783.99, 830.61, 880.00, 932.33, 987.77, 1046.50, 1108.73, 1174.66, 1244.51, 1318.51, 1396.91, 1479.98, 1567.98, 1661.22, 1760.00, 1864.66, 1975.53, 2093.00, 2217.46, 2349.32, 2489.02, 2637.02, 2793.83, 2959.96, 3135.96, 3322.44, 3520.00, 3729.31, 3951.07, 4186.01, 4434.92, 4698.64, 4978.03, 5274.04, 5587.65, 5919.91, 6271.93, 6644.88, 7040.00, };

unsigned int frequency = 44100;
unsigned int sampleFrequency = 0;
unsigned int audioBufferSize = 0;
unsigned int outputAudioBufferSize = 0;

const unsigned short sampleRange = 32767;

int gate;

struct Envelope {
    float a, d, s, r;
    float t;
    int _gate;
    char state; 
    void trigger() {
        _gate = 1;
        state = 0; 
        t = 0;
    }
    void process(float * samples, int count, int gate) {
        float m = 0;
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
                    m = t;
                    t += 1 / a / (float)frequency;
                    break;
                }
                t = 0;
                state++;
            case 1:
                if (d && state < 1) {
                    m = (1 - t * (1 - s));
                    t += 1 / d / (float)frequency;
                    if (t >= 0) {
                        state++;
                        t = 0;
                    }
                }
                t = 0;
                state++;
            case 2:
                m = s;
                break;
            case 3:
                if (r && t < 1) {
                    m = s * (1-t);
                    t += 1 / r / frequency;
                    break;
                }
                t = 0;
                state = 4;
            case 4:
                m = 0;
            }
            samples[i] *= m;
        }
    }
};

const float fade = 0.05;

struct Frequency {
    float freq;
    float freq0;
    float glide;
    float t;
    void generate(float * samples, int count) {
        if (freq != freq0 && t >= 1) {
            t = 0;
        }
        for (int i = 0; i < count; i++) {
            float f = (freq * t + freq0 * (1 - t));
            t += 1 / (float)frequency / glide;
            samples[i] = f;
        }
    }
};

struct Freq2Signal {
    float sine;
    float tri;
    float square;
    float t;
    void generate(float * samples, int count) {
        for (int i=0; i<count; i++) {
            float h = sin(2.0 * M_PI * t) * sine;
            h += (2.0 * t - 1.0) * tri;
            h += (t < 0.5 ? -1 : 1) * square;
            t += samples[i] / frequency;
            samples[i] = h;
            t = fmod(t, 1);
        }
    }
};

struct Signal {
    float freq;
    float sine;
    float tri;
    float square;
    float freq0;
    float t;
    float fadet;
    void generate(float * samples, int count) {
        if (freq != freq0 && fadet>=1) {
            fadet = 0;
        }
        for (int i=0; i<count; i++) {
            float f = (freq * fadet + freq0 * (1 - fadet));
            float h = sin(2.0 * M_PI * t) * sine;
            h += (2.0 * t - 1.0) * tri;
            h += (t < 0.5 ? -1 : 1) * square;
            t += f / frequency;
            fadet += 1 / (float)frequency / fade;
            if (fadet >= 1) {
                fadet = 1;
                freq0 = freq;
            }

            samples[i] = h;
            t = fmod(t, 1);
        }
    }
};

Frequency freqGen{440, 440, 0.1};
Freq2Signal freq2Sig{ 1 };
Signal voice {440, 1, 0, 0};
Envelope env{0.1, 0.1, 0.8, 1.0, 0, 0, 4};

void fillAudio(void *unused, Uint8 *stream, int len) {
    voice.generate((float*)stream, len/4);
    env.process((float*)stream, len/4, gate);
}

void noteOn(SDL_Event & e, int note) {
    if (e.key.repeat) return;
    gate++;
    //voice.freq = notes[note];

    env.trigger();
}

int main(int argc, char *argv[]) {
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

    SDL_PauseAudio(0);

    bool running = true;
    SDL_Event event;
    while (running) {
        while (SDL_PollEvent(&event)) {
            switch (event.type) {
            case SDL_KEYDOWN:
                switch (event.key.keysym.sym) {
                case SDLK_ESCAPE: running = false; break;
                case SDLK_a: noteOn(event, noteA4); break;
                case SDLK_s: noteOn(event, noteA4+2); break;
                case SDLK_d: noteOn(event, noteA4+4); break;
                case SDLK_f: noteOn(event, noteA4+5); break;
                case SDLK_g: noteOn(event, noteA4+7); break;
                case SDLK_h: noteOn(event, noteA4+9); break;
                case SDLK_j: noteOn(event, noteA4+11); break;
                case SDLK_k: noteOn(event, noteA4+12); break;
                case SDLK_F4: if (event.key.keysym.mod & KMOD_ALT) running = false; break;
                }
                break;
            case SDL_KEYUP:
                switch (event.key.keysym.sym) {
                case SDLK_a: gate--; break;
                case SDLK_s: gate--; break;
                case SDLK_d: gate--; break;
                case SDLK_f: gate--; break;
                case SDLK_g: gate--; break;
                case SDLK_h: gate--; break;
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
    }
    printf("quitting\n");
    SDL_Quit();
    return 0;
}
