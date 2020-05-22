#include <iostream>
#include <cmath>
#include <SDL2/SDL.h>
#include <SDL2/SDL_main.h>

/* linker options: -lmingw32 -lSDLmain -lSDL -mwindows */

using namespace std;

const int noteA0 = 0;
const int noteA4 = 48;
float notes[] = {27.50, 29.14, 30.87, 32.70, 34.65, 36.71, 38.89, 41.20, 43.65, 46.25, 49.00, 51.91, 55.00, 58.27, 61.74, 65.41, 69.30, 73.42, 77.78, 82.41, 87.31, 92.50, 98.00, 103.83, 110.00, 116.54, 123.47, 130.81, 138.59, 146.83, 155.56, 164.81, 174.61, 185.00, 196.00, 207.65, 220.00, 233.08, 246.94, 261.63, 277.18, 293.66, 311.13, 329.63, 349.23, 369.99, 392.00, 415.30, 440.00, 466.16, 493.88, 523.25, 554.37, 587.33, 622.25, 659.26, 698.46, 739.99, 783.99, 830.61, 880.00, 932.33, 987.77, 1046.50, 1108.73, 1174.66, 1244.51, 1318.51, 1396.91, 1479.98, 1567.98, 1661.22, 1760.00, 1864.66, 1975.53, 2093.00, 2217.46, 2349.32, 2489.02, 2637.02, 2793.83, 2959.96, 3135.96, 3322.44, 3520.00, 3729.31, 3951.07, 4186.01, 4434.92, 4698.64, 4978.03, 5274.04, 5587.65, 5919.91, 6271.93, 6644.88, 7040.00, };

unsigned int sampleFrequency = 0;
unsigned int audioBufferSize = 0;
unsigned int outputAudioBufferSize = 0;

const unsigned short sampleRange = 32767;

int gate;

struct Envelope {
    int a, d, r;
    float s;
    int t;
    int _gate;
    void process(Uint16 * samples, int count) {
        if (!_gate && gate) { t = -a; } // start attack
        if (!gate && gate) { t = 0; } // start release
        _gate = gate;
        for (int i=0; i<count; i++) {
            if (t < 0) { // attack
                samples[i] *= ((float)(a+t) / (float)t);
            } else if (!_gate && t > d) { // release
                samples[i] *= (s * (float)(t-d) / (float)r );
            } else if (gate && t == 0) { // sustain
                samples[i] *= s;
            } else if (gate && t > 0) { // decay
                float delta = (float)d - (float)t;
                float invs = 1.0 - s;
                samples[i] *= (1.0f - invs * delta);
            }
            if ((_gate && t < d) || !_gate) t = t + 1;
        }
    }
};

struct SineVoice {
    float freq;
    float gain;
    int _gate;
    int sample;
    float freq0;
    void generate(Uint16 * samples, int count) {
        if (freq0 != freq) { sample = 0; freq0 = freq; }
        if (_gate && !gate) { sample = 0; }
        _gate = gate;
        int period  = sampleFrequency / freq;
        for (int i=0; i<count; i++) {
            if (!_gate) { samples[i] = 0; continue; }
            float h = sin(2.0 * M_PI * (float)sample / (float)period);
            samples[i] = h * sampleRange * gain;
            sample = (sample + 1) % period;
        }
    }
};

SineVoice voice {440, 0.05};

void fillAudio(void *unused, Uint8 *stream, int len) {
    voice.generate((Uint16*)stream, len/2);
}

int main(int argc, char *argv[]) {
    if( SDL_Init(SDL_INIT_TIMER | SDL_INIT_AUDIO ) <0 ) {
        cout << "Unable to init SDL: " << SDL_GetError() << endl;
        return 1;
    }

    SDL_AudioSpec desired, obtained;
    desired.freq = 44100;
    desired.format = AUDIO_S16;
    desired.samples = 4096;
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
    if (obtained.format != AUDIO_S16) {
        fprintf(stderr, "did not get desired format of AUDIO_U16");
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
                case SDLK_a: gate++; voice.freq = notes[noteA4]; break;
                case SDLK_s: gate++; voice.freq = notes[noteA4+2]; break;
                case SDLK_d: gate++; voice.freq = notes[noteA4+4]; break;
                case SDLK_f: gate++; voice.freq = notes[noteA4+5]; break;
                case SDLK_g: gate++; voice.freq = notes[noteA4+7]; break;
                case SDLK_h: gate++; voice.freq = notes[noteA4+9]; break;
                case SDLK_j: gate++; voice.freq = notes[noteA4+11]; break;
                case SDLK_k: gate++; voice.freq = notes[noteA4+12]; break;
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
        SDL_Delay(1);
    }
    printf("quitting\n");
    SDL_Quit();
    return 0;
}
