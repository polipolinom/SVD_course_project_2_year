#include <gtest/gtest.h>

#include "graphic_tests.h"
#include "test_QR_decomposition.h"
#include "test_bidiagonalization.h"
#include "test_complex.h"
#include "test_matrix.h"
#include "test_orthonormalize.h"
#include "test_svd.h"
#include "test_vector.h"

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
