#include "gl/glew.h"
#include "glfw/glfw3.h"
#include <iostream>
#include <chrono>
#include <assert.h>
#include <vector>
#include <fstream>


const char *title = "test";

const char *vertCode = "\
#version 450\n\
layout(location = 0) in vec3 inPos;\n\
layout(location = 1) in vec3 inColor;\n\
layout(location = 0) out vec2 outTexcoord; \n\
out gl_PerVertex\n\
{\n\
    vec4 gl_Position;\n\
};\n\
void main() {\n\
  gl_Position = vec4(inPos, 1);\n\
  outTexcoord=2*vec2(gl_VertexID&1,gl_VertexID>>1);\n\
}\n\
";

const char *fragShaderPath = "shading.frag";

struct Vertex
{
    float pos[3];
    float color[3];
};

struct UBO
{
    float viewportSize[2];
    float mousePos[2]; //[0,1]
    int mouseButtonAction[3];
    float cameraDist{20};
} uboData{};

std::chrono::high_resolution_clock::time_point startTimes[2];

// triangle 1
float scale = 2.;
float offset = 1.;

std::vector<Vertex> vertices{
    {{-1 * scale + offset, -1 * scale + offset, 0}, {0, 1, 0}},
    {{1 * scale + offset, -1 * scale + offset, 0}, {0, 1, 0}},
    {{-1 * scale + offset, 1 * scale + offset, 0}, {0, 1, 0}},
};

GLuint createShader(GLenum stage, const char *code)
{
    GLuint ret = glCreateShader(stage);
    glShaderSource(ret, 1, &code, NULL);
    glCompileShader(ret);
    glReleaseShaderCompiler();

    GLint success;
    glGetShaderiv(ret, GL_COMPILE_STATUS, &success);
    if (!success)
    {
        GLint len;
        glGetShaderiv(ret, GL_INFO_LOG_LENGTH, &len);
        std::string log;
        log.resize(static_cast<size_t>(len));
        glGetShaderInfoLog(ret, len, NULL, &log[0]);
        std::cout << log << std::endl;
        glDeleteShader(ret);
        assert(0 && "failed to compile shader");
    }
    return ret;
}

GLuint createProgram(GLuint *shaders, size_t count)
{
    GLuint ret = glCreateProgram();
    for (size_t i = 0; i < count; ++i)
    {
        glAttachShader(ret, shaders[i]);
    }
    glLinkProgram(ret);
    return ret;
}

static void key_callback(GLFWwindow *window, int key, int scancode, int action, int mods);
static void mouse_button_callback(GLFWwindow *window, int button, int action, int mods);
static void scroll_callback(GLFWwindow *window, double xoffset, double yoffset);
static void framebuffer_size_callback(GLFWwindow *window, int width, int height);
static void cursor_position_callback(GLFWwindow *window, double xpos, double ypos);

GLFWwindow *window;
GLuint vbo;
GLuint program;
GLuint vao;
GLuint ubo;
bool updateUboFlag{true};

void createWindow()
{
    glfwInit();
    glfwWindowHint(GLFW_RESIZABLE, GL_TRUE);
    window = glfwCreateWindow(800, 600, title, nullptr, nullptr);
    glfwSetKeyCallback(window, key_callback);
    glfwSetMouseButtonCallback(window,mouse_button_callback);
    glfwSetScrollCallback(window, scroll_callback);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    glfwSetCursorPosCallback(window, cursor_position_callback);
    int width,height;
    glfwGetFramebufferSize(window, &width, &height);
    uboData.viewportSize[0] = width;
    uboData.viewportSize[1] = height;
    glfwMakeContextCurrent(window);
}
void configOpengl()
{
    GLenum err = glewInit();
    if (GLEW_OK != err)
    {
        printf("Error: %s\n", glewGetErrorString(err));
        assert(0 && "failed to init glew");
    }
    std::cout << glGetString(GL_VENDOR) << std::endl;
    std::cout << glGetString(GL_RENDERER) << std::endl;
    std::cout << glGetString(GL_VERSION) << std::endl;
    std::cout << glGetString(GL_SHADING_LANGUAGE_VERSION) << std::endl;
    glClearColor(0, 0.02, 0, 1);
    glEnable(GL_FRAMEBUFFER_SRGB);
}

void createVBO()
{
    glGenBuffers(1, &vbo);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, sizeof(Vertex) * vertices.size(), vertices.data(), GL_DYNAMIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}
void createProgram()
{
    std::fstream fs(fragShaderPath);
    if (!fs)
    {
        assert(0 && "failed to load shader");
    }
    fs.seekg(0, std::ios::end);
    auto len = fs.tellg();
    fs.seekg(0, std::ios::beg);
    std::string code;
    code.resize(len);
    fs.read(&code[0], len);
    fs.close();

    auto vertshader = createShader(GL_VERTEX_SHADER, vertCode);
    auto fragshader = createShader(GL_FRAGMENT_SHADER, code.data());
    GLuint shades[]{vertshader, fragshader};
    program = createProgram(shades, 2);
}
void createVAO()
{
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glEnableVertexAttribArray(0);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), 0);
    auto offset = offsetof(Vertex, color);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_TRUE, 6 * sizeof(float), (void *)offset);
    glBindVertexArray(0);
}
void createUBO()
{
    glGenBuffers(1, &ubo);
    glBindBuffer(GL_ARRAY_BUFFER, ubo);
    glBufferData(GL_ARRAY_BUFFER, sizeof(UBO), &uboData, GL_DYNAMIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void updateUBO()
{
    glBindBuffer(GL_UNIFORM_BUFFER, ubo);
    void *data = glMapBuffer(GL_UNIFORM_BUFFER, GL_WRITE_ONLY);
    memcpy(data, &uboData, sizeof(UBO));
    glUnmapBuffer(GL_UNIFORM_BUFFER);
}

void draw()
{
    if (updateUboFlag)
    {
        updateUBO();
        updateUboFlag = false;
    }
    
    glBindBufferRange(GL_UNIFORM_BUFFER, 0, ubo, 0, sizeof(UBO));
    glClear(GL_COLOR_BUFFER_BIT);
    glUseProgram(program);
    glBindVertexArray(vao);
    glDrawArrays(GL_TRIANGLES, 0, vertices.size());
}

void mainloop()
{
    while (!glfwWindowShouldClose(window))
    {
        glfwPollEvents();
        draw();
        glfwSwapBuffers(window);
    }
}
void key_callback(GLFWwindow *window, int key, int scancode, int action, int mods)
{
}
void mouse_button_callback(GLFWwindow *window, int button, int action, int mods)
{
    if (button == GLFW_MOUSE_BUTTON_LEFT)
        uboData.mouseButtonAction[0] = (action == GLFW_PRESS);
    if (button == GLFW_MOUSE_BUTTON_MIDDLE)
        uboData.mouseButtonAction[1] = (action == GLFW_PRESS);
    if (button == GLFW_MOUSE_BUTTON_RIGHT)
        uboData.mouseButtonAction[2] = (action == GLFW_PRESS);

    updateUboFlag = true;
}
void scroll_callback(GLFWwindow *window, double xoffset, double yoffset)
{
    uboData.cameraDist += yoffset;
    updateUboFlag = true;
}
void framebuffer_size_callback(GLFWwindow *window, int width, int height)
{
    while (width == 0 || height == 0)
    {
        glfwWaitEvents();
        glfwGetFramebufferSize(window, &width, &height);
    }
    uboData.viewportSize[0] = width;
    uboData.viewportSize[1] = height;
    updateUboFlag = true;
    glViewport(0, 0, width, height);
}
void cursor_position_callback(GLFWwindow *window, double xpos, double ypos)
{
    if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS)
    {
        uboData.mousePos[0] = xpos / uboData.viewportSize[0];
        uboData.mousePos[1] = 1 - ypos / uboData.viewportSize[1];
        updateUboFlag = true;
    }
}

int main()
{
    createWindow();
    configOpengl();
    createProgram();
    createVBO();
    createVAO();
    createUBO();
    mainloop();
}
