#pragma once

#include <sys/time.h>

class timer {
private:
    struct timeval m_StartTime, m_EndTime;

public:
    timer() {
        this->start();
    }

    void start() {
        gettimeofday(&m_StartTime, nullptr);
    }

    long millis() {
        long useconds, seconds, mseconds;
        gettimeofday(&m_EndTime, nullptr);
        useconds = m_EndTime.tv_usec - m_StartTime.tv_usec;
        seconds = m_EndTime.tv_sec - m_StartTime.tv_sec;
        mseconds = ((seconds) * 1000 + useconds/1000.0 + 0.5);
        return mseconds;
    }

    long seconds() {
        return millis() / 1000;
    }
};