#include <stdio.h>
#ifdef _WIN32
#include <windows.h>
#endif

#include "Timer.h"
#define USE_PERFORMANCE_COUNTER

#ifdef USE_PERFORMANCE_COUNTER

double TimerClass::SecondsPerTick = 0.0;

TimerClass::TimerClass()
{
	m_start = 0;
	m_end = 0;

	if (SecondsPerTick == 0.0){
		_LARGE_INTEGER now;
		QueryPerformanceFrequency (&now);
		SecondsPerTick = 1.0 / (double) now.QuadPart;
  } //end-if
}


void TimerClass::StartTiming (void)
{
	_LARGE_INTEGER now;
	if (QueryPerformanceCounter(&now) != 0)
	{
		m_start = now.QuadPart;
	}
	else
	{
		m_start = (TIMER_TimeT) GetTickCount(); // Fallback
	}
}

void TimerClass::StopTiming (void)
{
	_LARGE_INTEGER now;
	if (QueryPerformanceCounter(&now) != 0)
	{
		m_end = now.QuadPart;
	}
	else
	{
		m_end = (TIMER_TimeT) GetTickCount(); // Fallback
	}
}


// Returns the current time in miliseconds
__int64 TimerClass::GetTime(){
  static _LARGE_INTEGER now;
	QueryPerformanceCounter(&now);

	if (SecondsPerTick == 0.0) SecondsPerTick = 1.0 / (double) now.QuadPart;   // Initialize seconds per tick.
  static double MiliSecondsPerTick = SecondsPerTick*1e3;
  return (__int64)(MiliSecondsPerTick * now.QuadPart);
} //end-GetTime


double TimerClass::SecondsBetween ()
{
	return SecondsPerTick * double (m_end - m_start);
}

#elif defined (_WIN32)
double TimerClass::SecondsPerTick = 0.001;		// 1/1000

void TimerClass::StartTiming (void)
{
	m_start = (TIMER_TimeT) GetTickCount();
}

void TimerClass::PauseTiming(void) {
	//m_start = (TIMER_TimeT) GetTickCount();
}

void TimerClass::ContinueTiming(void) {
	//m_start = (TIMER_TimeT) GetTickCount();
}

void TimerClass::StopTiming (void)
{
	m_end = (TIMER_TimeT) GetTickCount();
}

double TimerClass::SecondsBetween ()
{
	return 0.001 * double (m_end - m_start);
}

#else // UNIX

TimerClass::TimerClass()
{
	m_start.tv_sec  = 0;
	m_start.tv_usec = 0;
	m_end.tv_sec    = 0;
	m_end.tv_usec   = 0;
}

void TimerClass::StartTiming (void)
{
	gettimeofday (&m_start, NULL);
}

void TimerClass::StopTiming (void)
{
	gettimeofday (&m_end, NULL);
}

double TimerClass::SecondsBetween (void)
{
	double t;

	t = double (m_end.tv_usec - m_start.tv_usec) * 1.e-6;
	t += double (m_end.tv_sec - m_start.tv_sec);

	return t;
}
#endif


void TimerClass::OutputTime (const char * where)
{
//	double t = SecondsBetween ();
//	char buffer[1000];

//	sprintf (buffer, "TIMER: %8.3f secs for %s\n", t, where);
//	OutputDebugString (buffer);
}

void TimerClass::OutputSpeed (const char * where, char * units, double mult)
{
	double t = SecondsBetween ();
//	char buffer[1000];

//	if (t != 0.0) t = mult / t;

//	sprintf (buffer, "TIMER: %8.3f %s for %s\n", t, (units ? units : "per second"), where);
//	OutputDebugString (buffer);
}
