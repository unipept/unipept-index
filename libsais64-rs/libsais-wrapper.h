#include "libsais-packed/libsais/include/libsais16x64.h"
#include "libsais-packed/libsais/include/libsais32x64.h"
#include "libsais-packed/libsais/include/libsais64.h"


int64_t libsais16x64(const uint16_t * T, int64_t * SA, int64_t n, int64_t fs, int64_t * freq);

int64_t libsais32x64(const uint32_t * T, int64_t * SA, int64_t n, int64_t k, int64_t fs, int64_t * freq);

int64_t libsais64(const uint8_t * T, int64_t * SA, int64_t n, int64_t fs, int64_t * freq);