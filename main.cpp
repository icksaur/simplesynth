#include <cmath>
#include <ctime>
#include <SDL2/SDL.h>
#include <SDL2/SDL_main.h>

#include "GLObject.h"
#include "glxw.h"
#include "gui.h"
#include "synth.h"

using namespace std;

const int windowResolution = 512;
int gate;

int graph[windowResolution];
bool graphDirty = false;
void renderWaveform(float * samples, int count) {
    const float min = -1.0f;
    const float max = 1.0f;
    const float range = max - min;
    if (!graphDirty) return;
    graphDirty = false;
    int width = count < windowResolution ? count : windowResolution;
    for (int i = 0; i < width; i++) {
        graph[i] = windowResolution - (int)((samples[i] - min) / range * (float)windowResolution);
    }
}

Note noteGen{ 0.001f, 4, 4 };
Frequency noteFreq;
ModulatingSignal modsig { 0, 1, 0, 0 };
FixedSignal fixedsig { 0, 1, 0, 0 };
Envelope env{0.05f, 0, 1, 0.05f, 0, 0, 4};
LowpassFilter lowpassFilter;
AllPassFilter allPassFilter{ 0.5f };
AllPassDelay allPassDelay{ -0.5f, 22050 };
Delay delay{ 0.5f, 0.5f };
ExplicitDelay ed1{ 229, 0.5f };
ExplicitDelay ed2{ 1069, 0.5f };
ExplicitDelay ed3{ 3181, 0.5f };
ExplicitDelay ed4{ 6053, 0.5f };
ExplicitDelay ed5{ 7919, 0.5f };
StereoReverb stereoReverb;
MoogFilter moogFilter(0.5f, 0.0f);
VocalFilter vocalFilter;
LFO lfo{ 0.02f, 0.5f };
Fader fader;

Sampler sampler(frequency * 10);
SamplerWriteHead writeHead(sampler);
SamplerReadHead readHead(sampler);

int mode = 0;
bool isPlaying = false;
bool isWriting = false;

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
            fader.process(sample, stereoSamples[i]);
        } else if (mode == 1) {
            vocalFilter.process(sample);
            stereoReverb.process(sample, stereoSamples[i]);
            lowpassFilter.process(stereoSamples[i].r);
            lowpassFilter.process(stereoSamples[i].l);
        } else {
            fader.process(sample, stereoSamples[i]);
        }

        if (isWriting) writeHead.write((float*)(&stereoSamples[i]), 2);
        if (isPlaying && !isWriting) readHead.process((float*)(&stereoSamples[i]), 2);
    }
    renderWaveform((float*)stream, len);
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
    srand((unsigned)time(0)); // for noise generator
    int xwindowResolution = windowResolution, ywindowResolution = windowResolution;
    SDL_Init(SDL_INIT_VIDEO | SDL_INIT_AUDIO);
    SDL_Window * window = SDL_CreateWindow("test", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, xwindowResolution, xwindowResolution, SDL_WINDOW_OPENGL | SDL_WINDOW_ALLOW_HIGHDPI);
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

    unsigned graphBytesSize = windowResolution * windowResolution * 4;
    char * graphBytes = new char[graphBytesSize];
    unsigned * graphTexels = (unsigned*)graphBytes;
    memset(graphBytes, 0, graphBytesSize);
    ArrayTexture graphTexture(windowResolution, 1, graphBytes, graphBytesSize);

    GUI gui;
    int octave = noteA4;
    bool running = true;

    SDL_Event event;
    while (running) {
        while (SDL_PollEvent(&event)) {
            switch (event.type) {
            case SDL_KEYDOWN:
                if (event.key.repeat == 1) break;
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
                case SDLK_r: writeHead.writeHead = 0; isWriting = true; break;
                case SDLK_p: readHead.readHead = 0; isPlaying = true; break;
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
                case SDLK_r: isWriting = false; break;
                case SDLK_p: readHead.readHead = 0; isPlaying = false; break;
                }
                break;
            case SDL_QUIT:
                running = false;
                break;
            case SDL_MOUSEMOTION: gui.pushMouseMotion(event.motion.x, windowResolution - event.motion.y); break;
            case SDL_MOUSEBUTTONDOWN: if (event.button.button == SDL_BUTTON_LEFT) gui.pushMouseDown(0, event.button.x, windowResolution - event.button.y); break;
            case SDL_MOUSEBUTTONUP: if (event.button.button == SDL_BUTTON_LEFT) gui.pushMouseUp(0, event.button.x, windowResolution - event.button.y); break;
            }
        }

        memset(graphTexels, 0, sizeof(unsigned) * windowResolution * windowResolution);
        for (int i = 0; i < windowResolution; i++) {
            for (int j = 0; j < windowResolution; j++) {
                if (graph[i] < j) {
                    graphTexels[windowResolution * j + i] = 0xFFFFFFFF;
                }
            }
        }
        glActiveTexture(GL_TEXTURE0 + 0);
        graphTexture.bind();
        glTexSubImage3D(GL_TEXTURE_2D_ARRAY, 0, 0, 0, 0, windowResolution, windowResolution, 1, GL_RGBA, GL_UNSIGNED_BYTE, graphTexels);

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
        static float formants = 1.0f;

        const char * filterNames[3]{ "low", "high", "band" };

        switch (mode) {
        case 0:
            if (param(gui, formants, "formants")) {
                vocalFilter.formants = (int)((formants * 4.0f) + 1.0f);
            }
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
        if (isWriting) gui.label(0, elementY+step/2, "rec");
        else if (isPlaying) gui.label(0, elementY+step/2, "play");
        gui.end();
        gui.render();
        graphDirty = true;
        SDL_GL_SwapWindow(window);
        SDL_Delay(1000/60);
    }

    return 0;
}
