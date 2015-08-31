#ifndef GETTIMEMS64_HPP_INCLUDED
#define GETTIMEMS64_HPP_INCLUDED

#include <Snap.h>

#ifdef _WIN32
#include <Windows.h>
#else
#include <sys/time.h>
#include <ctime>
#endif

/* Remove if already defined */
//typedef long long int64; typedef unsigned long long uint64;

/* Returns the amount of milliseconds elapsed since the UNIX epoch. Works on both
 * windows and linux. */

uint64 GetTimeMs64();

#endif // GETTIMEMS64_HPP_INCLUDED