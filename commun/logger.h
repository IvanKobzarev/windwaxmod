#pragma once

#include <fstream>
#include <iostream>
#include <ctime>
#include <stdarg.h>

using namespace std;

class Logger {
  public:
    Logger(char* filename);
    ~Logger();
    void log(char* logline);
    void logp(char* logline);
    void logprintf(char *str,...);
    void close();
   private:
    ofstream m_stream;
};

void LoggerInit(char* filename);

void LayoutLoggerInit(char* filename);

Logger* getLogger();

Logger* getLayoutLogger();

char* getLogFileName();

