#ifndef TIMER_C_H
#define TIMER_C_H

#ifdef __cplusplus
extern "C" {
#endif

void swTimerStart(const char* name);
void swTimerEnd(const char* name);
void swTimerPrintAll();
void swTimerPrint(const char* name);
void swTimerMaxSum();

//- interface for Fortran
void swtimerstart_(const char* name);
void swtimerend_(const char* name);
void swtimerprintall_();
void swtimerprint_(const char* name);
void swtimermaxsum_();

#ifdef __cplusplus
}
#endif

#endif
