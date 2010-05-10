#include "logger.h"

Logger* logger=NULL, *layoutLogger=NULL;

void LoggerInit(char* filename) {
	logger = new Logger(filename);
}

void LayoutLoggerInit(char* filename) {
	layoutLogger = new Logger(filename);
        printf ("\n Log to %s", filename);
}


Logger* getLogger() {
	return logger;
}

Logger* getLayoutLogger() {
	return layoutLogger;
}


Logger::Logger(char* filename) {
  m_stream.open(filename);
}

void Logger::log(char* logline){
  m_stream << logline << endl;
}

void Logger::logp(char* logline){
  m_stream << logline << endl;
  cout << logline << endl;
}


void Logger::logprintf(char *str,...)
{
    char logline[1000];
    va_list arglist;
    va_start(arglist, str);
    vsprintf(logline, str, arglist);
    va_end(arglist);
    m_stream << logline;// << endl;
    cout << logline;;// << endl;
}


Logger::~Logger(){
  m_stream.close();
}

void Logger::close(){
  m_stream.close();
}


char* getLogFileName() {
    char filename[255];
    char date[10], time[10];
    _strdate(date);
    //printf ("\nDATE=[%s]", date);
    _strtime(time);
    //printf ("\nTIME=[%s]", time);
    for (int i = 0; i < 8; i++) if (date[i] == '/') date[i]='_';
    for (int i = 0; i < 8; i++) if (time[i] == ':') time[i]='_';
    //printf ("\n.DATE=[%s]", date);
    //printf ("\n.TIME=[%s]", time);
    sprintf (filename, "log%s-%s.txt",  date, time);
    //printf ("\n filename=%s", filename);
    return filename;
}