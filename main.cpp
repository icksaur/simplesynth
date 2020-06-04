#include <iostream>
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

SDL_Renderer * renderer;
int graph[512];
void renderState(float * samples, int count, float min = -1.0f, float max = 1.0f) {
    int width = count < 512 ? count : 512;
    float range = max - min;
    for (int i = 0; i < width; i++) {
        graph[i] = 512 - (samples[i] - min) / range * 512;
    }
}

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
    float k;
    float r;
    void process(float * stream, int count) {
        for (int i = 0; i < count; i++) {
            stream[i] = (stream[i] - r) * k + r;
            r = stream[i];
        }
    }
};

const int delayBufferSize = 88211; // prime larger than 2 * 44100

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
    void process(float * samples, int count) {
        delay = delay > (frequency - 1) ? (frequency - 1) : delay;
        for (int i = 0; i < count; i++) {
            process(samples[i]);
        }
    }
};

struct AllPassDelay {
    float k;
    int delay;
    float buffer[frequency];
    int read;
    void process(float * stream, int count) {
        for (int i = 0; i < count; i++) {
            float r1 = buffer[read] + stream[i] * k; // peek in delay
            stream[i] = r1 * -k + buffer[read];
            buffer[(read + delay) % frequency] = r1 * 0.5;
            ++read %= frequency;
        }
    }

};

// Not sure what this does. It's a building block for other things.
struct AllPassFilter {
    float k;
    float r;
    void process(float * stream, int count) {
        for (int i = 0; i < count; i++) {
            float r1 = stream[i] + r * k;
            stream[i] = r1 * -k + r;
            r = r1;
        }
    }
};

struct WindowFilter {
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

int baseDelay = 441;

struct Reverb {
    ExplicitDelay delays[8];
    ExplicitDelay tapDelays[4];
    float accum[1024];
    float output[1024];
    float krt;
    int phase;
    Reverb() {
        delays[0].delay = baseDelay + 29; delays[0].decay = 0.3;
        delays[1].delay = baseDelay + 601; delays[1].decay = 0.5;
        delays[2].delay = baseDelay + 1291; delays[2].decay = 0.7;
        delays[3].delay = baseDelay + 2053; delays[3].decay = 0.9;
        delays[4].delay = baseDelay + 2819; delays[4].decay = 0.9;
        delays[5].delay = baseDelay + 3643; delays[5].decay = 0.7;
        delays[6].delay = baseDelay + 4493; delays[6].decay = 0.5;
        delays[7].delay = baseDelay + 5387; delays[7].decay = 0.3;

        tapDelays[0].delay = baseDelay + 29; tapDelays[0].decay = 0.1;
        tapDelays[1].delay = baseDelay + 601; tapDelays[1].decay = 0.1;
        tapDelays[2].delay = baseDelay + 1291; tapDelays[2].decay = 0.1;
        tapDelays[3].delay = baseDelay + 2053; tapDelays[3].decay = 0.1;

        memset(accum, 0, 1024 * sizeof(float));
        phase = 0;
        krt = 0.1;
    }

    void process(float * samples, int count) {
        memset(output, 0, 1024 * sizeof(float));
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < count; j++) accum[j] += samples[j];
            delays[2*i+0].process(accum, count);
            delays[2*i+1].process(accum, count);
            for (int j = 0; j < count; j++) output[j] += accum[j];
            for (int j = 0; j < count; j++) accum[j] *= krt;
        }
        for (int i = 0; i < count; i++) samples[i] = output[i];
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
WindowFilter windowFilter;
BasicFilter basicFilter{ 0.5 };
AllPassFilter allPassFilter{ 0.5 };
AllPassDelay allPassDelay{ -0.5, 22050 };
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
    //combinedReverb.process(FLOATSTREAM);
    //delay.process(FLOATSTREAM);
    if (mode == 0) {
        reverb.process(FLOATSTREAM); // delay.process(FLOATSTREAM);
        basicFilter.process(FLOATSTREAM);
    } else if (mode == 1)
        allPassDelay.process(FLOATSTREAM);
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

void prepVAO() {
    ShaderProgram program(fullScreenTexturedQuadSource);

    VAOBuilder builder(1);
    builder.addFloats(0, 0, 2).addFloats(0, 1, 2);
    VAO vao(program, builder);

    const float w = 0.5f;
    float vertices[] = { -w,-w, 0,0, -w,w, 0,1, w,w, 1,1, w,-w, 1,0 };
    vao.vertexData(vertices, sizeof(vertices));
}

int main(int argc, char *argv[]) {
    srand(time(0));
    int xrez = 512, yrez = 512;
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
    desired.channels = 1;
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

    unsigned graphBytesSize = 512 * 512 * 4;
    char * graphBytes = new char[graphBytesSize];
    unsigned * graphTexels = (unsigned*)graphBytes;
    memset(graphBytes, 0, graphBytesSize);
    ArrayTexture graphTexture(512, 1, graphBytes, graphBytesSize);

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
            case SDL_MOUSEMOTION: gui.pushMouseMotion(event.motion.x, 512 - event.motion.y); break;
            case SDL_MOUSEBUTTONDOWN: if (event.button.button == SDL_BUTTON_LEFT) gui.pushMouseDown(0, event.button.x, 512-event.button.y); break;
            case SDL_MOUSEBUTTONUP: if (event.button.button == SDL_BUTTON_LEFT) gui.pushMouseUp(0, event.button.x, 512-event.button.y); break;
            }
        }

        for (int i = 0; i < 512; i++) {
            for (int j = 0; j < 512; j++) {
                unsigned texel = graph[i] < j ? 0xFFFFFFFF : 0x00000000;
                graphTexels[512 * j + i] = texel;
            }
        }
        glActiveTexture(GL_TEXTURE0 + 0);
        graphTexture.bind();
        glTexSubImage3D(GL_TEXTURE_2D_ARRAY, 0, 0, 0, 0, 512, 512, 1, GL_RGBA, GL_UNSIGNED_BYTE, graphTexels);

        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        program.use();
        vao.bind();
        glDrawArrays(GL_TRIANGLE_FAN, 0, 4);

        static int x, y;
        static char krtString[12];
        gui.start(x, y, 300, 200, "params");
        switch (mode) {
        case 0:
            gui.slider(reverb.krt, 0, 0, 150, 20);
            sprintf_s(krtString, "%.2f krt", reverb.krt);
            gui.label(150, 10, krtString);
            gui.slider(basicFilter.k, 0, 25, 100, 20);
            gui.label(150, 35, "filter");
            if (gui.button(0, 50, 100, 20, "panic")) {
                memset(reverb.accum, 0, sizeof(float) * 1024);
                for (int i = 0; i < 8; i++) {
                    memset(reverb.delays[i].buffer, 0, sizeof(float) * 88211);
                }
            }
        }
        gui.end();
        gui.render();
        SDL_GL_SwapWindow(window);
        SDL_Delay(1000 / 120);
    }

exit:
    SDL_GL_DeleteContext(context);
    SDL_DestroyWindow(window);
    return 0;
}
