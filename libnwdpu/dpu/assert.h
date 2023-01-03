/*
 * Copyright 2022 - UPMEM
 */

#ifndef D78E8DA6_4F39_439C_BBFA_22ECA8BA9819
#define D78E8DA6_4F39_439C_BBFA_22ECA8BA9819

#define assert(expr) \
    if (!(expr))     \
    __asm__ volatile("fault 3")

#endif /* D78E8DA6_4F39_439C_BBFA_22ECA8BA9819 */
