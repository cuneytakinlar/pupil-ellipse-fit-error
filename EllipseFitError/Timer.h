#ifndef TIMER_HEADER_H
#define TIMER_HEADER_H 1

#ifdef _WIN32
#define USE_PERFORMANCE_COUNTER
#endif

#ifdef USE_PERFORMANCE_COUNTER

typedef __int64 TIMER_TimeT;
TIMER_TimeT TimerTimeStamp (void);

#elif define (_WIN32)

typedef long TIMER_TimeT;
inline TIMER_TimeT TimerTimeStamp (void)
{
	return (TIMER_TimeT) GetTickCount();
}
#else // Assume UNIX

#include <sys/time.h>

typedef struct timeval	TIMER_TimeT;

#endif


/** \class TimerClass
 **
 ** \brief A class for reading time and measuring elapsed time.
 **
 */


class TimerClass
{
private:
	TIMER_TimeT			m_start;
	TIMER_TimeT			m_end;
	static double		SecondsPerTick;

public:
  /// Constructor
	TimerClass ();
  /// Start timing.
	void StartTiming (void);
  /// Pause timing
  void PauseTiming(void);
  /// Continue timing
  void ContinueTiming(void);
  /// Stop timing
	void StopTiming (void);
  /// Elapsed time
	double SecondsBetween (void);    

  /// Returns the current time in miliseconds
  static __int64 GetTime();
    
  /// prints the time into the given string.
	void OutputTime (const char * string);
	void OutputSpeed (const char * where, char * units = 0, double mult = 1.0);
};

#ifdef qTiming
#define	TIMER_DECLARE_TIMER(name)			TimerClass name
#define TIMER_START_TIMING(name)			name.StartTiming()
#define TIMER_STOP_TIMING(name)				name.StopTiming()
#define TIMER_OUTPUT_TIMING(name,string)	name.OutputTime (string)
#define TIMER_OUTPUT_SPEED(name,s1, s2)		name.OutputSpeed (s1, s2)
#else
#define TIMER_DECLARE_TIMER(name)			/* nothing */
#define TIMER_START_TIMING(name)			/* nothing */
#define TIMER_STOP_TIMING(name)				/* nothing */
#define TIMER_OUTPUT_TIMING(name, string)	/* nothing */
#define TIMER_OUTPUT_SPEED(name,s1, s2)		/* nothing */
#endif


#endif