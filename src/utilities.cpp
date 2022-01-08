#include "../include/utilities.h"

void atg_scs::freeArray(double *&data) {
    delete[] data;
    data = nullptr;
}
