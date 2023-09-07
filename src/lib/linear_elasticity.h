#pragma once

#include <memory>

#include "laclib.h"

std::unique_ptr<Matrix> linear_elasticity_modulus(double young, double poisson, bool plane_stress) {
    auto dd = Matrix::make_new(4, 4);
    if (plane_stress) {
        double c = young / (1.0 - poisson * poisson);
        dd->set(0, 0, c);
        dd->set(0, 1, c * poisson);
        dd->set(1, 0, c * poisson);
        dd->set(1, 1, c);
        dd->set(3, 3, c * (1.0 - poisson)); // Mandel: multiply by 2, so 1/2 disappears
    } else {
        double c = young / ((1.0 + poisson) * (1.0 - 2.0 * poisson));
        dd->set(0, 0, c * (1.0 - poisson));
        dd->set(0, 1, c * poisson);
        dd->set(0, 2, c * poisson);
        dd->set(1, 0, c * poisson);
        dd->set(1, 1, c * (1.0 - poisson));
        dd->set(1, 2, c * poisson);
        dd->set(2, 0, c * poisson);
        dd->set(2, 1, c * poisson);
        dd->set(2, 2, c * (1.0 - poisson));
        dd->set(3, 3, c * (1.0 - 2.0 * poisson)); // Mandel: multiply by 2, so 1/2 disappears
    }
    return dd;
}
