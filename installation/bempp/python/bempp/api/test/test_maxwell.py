"""Validation tests for Maxwell problems."""

from unittest import TestCase
import bempp.api
import numpy as np

#pylint: disable=invalid-name
#pylint: disable=missing-docstring
#pylint: disable=unused-argument
#pylint: disable=too-many-locals


class TestMaxwell(TestCase):
    """Maxwell validation tests."""

    def test_efie_unit_sphere_rt_functions(self):
        """Exterior EFIE problem un the unit sphere with RT functions."""

        # This script solves the Maxwell equations in the region exterior
        # to a bounded
        # object, with Dirichlet boundary conditions given by the exact solution
        # (satisfying the Silver-Mueller radiation conditions)
        #
        #     \vec u(\vec x) = h_1^{(1)}(k r) \hat phi,
        #
        # where (r, theta, phi) are the radial, zenith angle and azimuthal
        # spherical
        # coordinates in the system anchored at the point
        # (0.1, 0.1, 0.1), h_1^{(1)}(r)
        # is the spherical Hankel function of the first kind and order
        # 1 and \hat phi is
        # the unit vector oriented along d(\vec x)/d\phi.

        k = 1
        source = 0.1

        grid = bempp.api.shapes.regular_sphere(4)
        parameters = bempp.api.common.global_parameters()
        parameters.assembly.boundary_operator_assembly_type = 'dense'

        def eval_dirichlet_data(point, normal, domain_index, result):
            x, y, z = point - source
            r = np.sqrt(x**2 + y**2 + z**2)
            kr = k * r
            h1kr = (-1j - kr) * np.exp(1j * kr) / (kr * kr)
            scale = h1kr / r
            field = [-y * scale, x * scale, 0.]
            result[:] = np.cross(field, normal)

        def eval_exact_neumann_data(point, normal, domain_index, result):
            x, y, z = point - source
            r = np.sqrt(x**2 + y**2 + z**2)
            kr = k * r
            h1kr = (-1j - kr) * np.exp(1j * kr) / (kr * kr)
            h1kr_deriv = ((1. + 1j - 1j * kr) * (1. + 1j + kr) *
                          np.exp(1j * kr) / (kr * kr * r))
            xy_factor = (h1kr - r * h1kr_deriv) / (r * r * r)
            curl = [x * z * xy_factor,
                    y * z * xy_factor,
                    ((x * x + y * y + 2 * z * z) * h1kr + r *
                     (x * x + y * y) * h1kr_deriv) /
                    (r * r * r)]
            result[:] = np.cross(curl, normal) / (1j * k)

#        def eval_exact_solution(point):
#            x, y, z = point - source
#            r = np.sqrt(x**2 + y**2 + z**2)
#            kr = k * r
#            h1kr = (-1j - kr) * np.exp(1j * kr) / (kr * kr)
#            scale = h1kr / r
#            return np.array([-y * scale, x * scale, 0.])

        rt_space = bempp.api.function_space(grid, "RT", 0)
        nc_space = bempp.api.function_space(grid, "NC", 0)

        efie = bempp.api.operators.boundary.maxwell.electric_field(
            rt_space, rt_space, nc_space, k, parameters=parameters,
            use_projection_spaces=False)
        mfie = bempp.api.operators.boundary.maxwell.magnetic_field(
            rt_space, rt_space, nc_space, k, parameters=parameters,
            use_projection_spaces=False)
        ident = bempp.api.operators.boundary.sparse.identity(
            rt_space, rt_space, nc_space, parameters=parameters)

        dirichlet_grid_fun = bempp.api.GridFunction(
            rt_space, fun=eval_dirichlet_data)
        rhs = -(.5 * ident + mfie) * dirichlet_grid_fun

        sol = bempp.api.linalg.lu(efie, rhs)

        exact_solution = bempp.api.GridFunction(
            rt_space, fun=eval_exact_neumann_data)
        rel_error = (sol - exact_solution).l2_norm() / exact_solution.l2_norm()
        self.assertTrue(
            rel_error < 2E-2,
            msg="Actual error: {0}. Expected error: 2E-2".format(rel_error))

    def test_efie_calderon_unit_sphere(self):
        """Exterior EFIE problem with Calderon preconditioning."""

        # This script solves the Maxwell equations in the region exterior
        # to a bounded
        # object, with Dirichlet boundary conditions given by the exact solution
        # (satisfying the Silver-Mueller radiation conditions)
        #
        #     \vec u(\vec x) = h_1^{(1)}(k r) \hat phi,
        #
        # where (r, theta, phi) are the radial, zenith angle and azimuthal
        # spherical coordinates in the system anchored at the point
        # (0.1, 0.1, 0.1), h_1^{(1)}(r)
        # is the spherical Hankel function of the first kind and order
        # 1 and \hat phi is
        # the unit vector oriented along d(\vec x)/d\phi.

        k = 1
        source = 0.1

        from bempp.api.operators.boundary.maxwell import \
            calderon_electric_field

        grid = bempp.api.shapes.regular_sphere(3)
        efie_squared, efie = calderon_electric_field(
            grid, k)
        rwg_space = efie_squared.domain
        snc_space = bempp.api.function_space(grid, "B-SNC", 0)
        bc_space = bempp.api.function_space(grid, "BC", 0)

        def eval_dirichlet_data(point, normal, domain_index, result):
            x, y, z = point - source
            r = np.sqrt(x**2 + y**2 + z**2)
            kr = k * r
            h1kr = (-1j - kr) * np.exp(1j * kr) / (kr * kr)
            scale = h1kr / r
            field = [-y * scale, x * scale, 0.]
            result[:] = np.cross(field, normal)

        def eval_exact_neumann_data(point, normal, domain_index, result):
            x, y, z = point - source
            r = np.sqrt(x**2 + y**2 + z**2)
            kr = k * r
            h1kr = (-1j - kr) * np.exp(1j * kr) / (kr * kr)
            h1kr_deriv = ((1. + 1j - 1j * kr) * (1. + 1j + kr) *
                          np.exp(1j * kr) / (kr * kr * r))
            xy_factor = (h1kr - r * h1kr_deriv) / (r * r * r)
            curl = [x * z * xy_factor,
                    y * z * xy_factor,
                    ((x * x + y * y + 2 * z * z) * h1kr +
                     r * (x * x + y * y) * h1kr_deriv) /
                    (r * r * r)]
            result[:] = np.cross(curl, normal) / (1j * k)

#        def eval_exact_solution(point):
#            x, y, z = point - source
#            r = np.sqrt(x**2 + y**2 + z**2)
#            kr = k * r
#            h1kr = (-1j - kr) * np.exp(1j * kr) / (kr * kr)
#            scale = h1kr / r
#            return np.array([-y * scale, x * scale, 0.])

        mfie = bempp.api.operators.boundary.maxwell.magnetic_field(
            rwg_space, bc_space, snc_space, k)
        ident = bempp.api.operators.boundary.sparse.identity(
            rwg_space, bc_space, snc_space)

        dirichlet_grid_fun = bempp.api.GridFunction(
            rwg_space, fun=eval_dirichlet_data)
        rhs = - efie * (.5 * ident + mfie) * dirichlet_grid_fun

        sol, _, residuals = bempp.api.linalg.gmres(
            efie_squared, rhs,
            use_strong_form=True, return_residuals=True)

        exact_solution = bempp.api.GridFunction(
            rwg_space, fun=eval_exact_neumann_data)
        rel_error = (sol - exact_solution).l2_norm() / exact_solution.l2_norm()
        self.assertTrue(
            rel_error < 5E-2,
            msg="Actual error: {0}. Expected error: 5E-2".format(rel_error))
        self.assertTrue(
            len(residuals) < 7,
            msg="Needed {0} iterations to ".format(len(residuals)) +
            "solve the system. Expected not more than 6 iterations.")

    def test_efie_unit_sphere_rwg_functions(self):
        """Exterior EFIE problem on the unit sphere with RWG functions."""

        # This script solves the Maxwell equations in the region exterior
        # to a bounded
        # object, with Dirichlet boundary conditions given by the exact solution
        # (satisfying the Silver-Mueller radiation conditions)
        #
        #     \vec u(\vec x) = h_1^{(1)}(k r) \hat phi,
        #
        # where (r, theta, phi) are the radial, zenith angle and
        #  azimuthal spherical
        # coordinates in the system anchored at the point
        #  (0.1, 0.1, 0.1), h_1^{(1)}(r)
        # is the spherical Hankel function of the first kind
        #  and order 1 and \hat phi is
        # the unit vector oriented along d(\vec x)/d\phi.

        k = 1
        source = 0.1

        grid = bempp.api.shapes.regular_sphere(4)
        parameters = bempp.api.common.global_parameters()
        parameters.assembly.boundary_operator_assembly_type = 'dense'

        def eval_dirichlet_data(point, normal, domain_index, result):
            """Evaluate the Dirichlet data."""
            x, y, z = point - source
            r = np.sqrt(x**2 + y**2 + z**2)
            kr = k * r
            h1kr = (-1j - kr) * np.exp(1j * kr) / (kr * kr)
            scale = h1kr / r
            field = [-y * scale, x * scale, 0.]
            result[:] = np.cross(field, normal)

        def eval_exact_neumann_data(point, normal, domain_index, result):
            """Evaluate the exact Neumann data."""
            x, y, z = point - source
            r = np.sqrt(x**2 + y**2 + z**2)
            kr = k * r
            h1kr = (-1j - kr) * np.exp(1j * kr) / (kr * kr)
            h1kr_deriv = ((1. + 1j - 1j * kr) * (1. + 1j + kr) *
                          np.exp(1j * kr) / (kr * kr * r))
            xy_factor = (h1kr - r * h1kr_deriv) / (r * r * r)
            curl = [x * z * xy_factor,
                    y * z * xy_factor,
                    ((x * x + y * y + 2 * z * z) * h1kr + r *
                     (x * x + y * y) * h1kr_deriv) /
                    (r * r * r)]
            result[:] = np.cross(curl, normal) / (1j * k)

#        def eval_exact_solution(point):
#            x, y, z = point - source
#            r = np.sqrt(x**2 + y**2 + z**2)
#            kr = k * r
#            h1kr = (-1j - kr) * np.exp(1j * kr) / (kr * kr)
#            scale = h1kr / r
#            return np.array([-y * scale, x * scale, 0.])

        rwg_space = bempp.api.function_space(grid, "RWG", 0)
        snc_space = bempp.api.function_space(grid, "SNC", 0)

        efie = bempp.api.operators.boundary.maxwell.electric_field(
            rwg_space, rwg_space, snc_space, k, parameters=parameters,
            use_projection_spaces=False)
        mfie = bempp.api.operators.boundary.maxwell.magnetic_field(
            rwg_space, rwg_space, snc_space, k, parameters=parameters,
            use_projection_spaces=False)
        ident = bempp.api.operators.boundary.sparse.identity(
            rwg_space, rwg_space, snc_space, parameters=parameters)

        dirichlet_grid_fun = bempp.api.GridFunction(
            rwg_space, fun=eval_dirichlet_data)
        rhs_coeffs = -(.5 * ident.weak_form() + mfie.weak_form()
                      ) * dirichlet_grid_fun.coefficients

        from scipy.linalg import solve
        sol_coefficients = solve(
            bempp.api.as_matrix(efie.weak_form()), rhs_coeffs)
        sol = bempp.api.GridFunction(rwg_space, coefficients=sol_coefficients)

        exact_solution = bempp.api.GridFunction(
            rwg_space, fun=eval_exact_neumann_data)
        rel_error = (sol - exact_solution).l2_norm() / exact_solution.l2_norm()
        self.assertTrue(
            rel_error < 2E-2,
            msg="Actual error: {0}. Expected error: 2E-2".format(rel_error))

    def test_maxwell_potential_operators(self):
        """Exterior EFIE problem on the unit sphere with RWG functions."""

        # This script solves the Maxwell equations in the region exterior to a
        # bounded
        # object, with Dirichlet boundary conditions given by the exact solution
        # (satisfying the Silver-Mueller radiation conditions)
        #
        #     \vec u(\vec x) = h_1^{(1)}(k r) \hat phi,
        #
        # where (r, theta, phi) are the radial, zenith angle and azimuthal
        # spherical
        # coordinates in the system anchored at the point
        # (0.1, 0.1, 0.1), h_1^{(1)}(r)
        # is the spherical Hankel function of the first kind and order
        # 1 and \hat phi is
        # the unit vector oriented along d(\vec x)/d\phi.

        k = 1
        source = 0.1

        grid = bempp.api.shapes.regular_sphere(5)
        parameters = bempp.api.common.global_parameters()
        parameters.assembly.boundary_operator_assembly_type = 'dense'
        parameters.assembly.potential_operator_assembly_type = 'dense'

        def eval_dirichlet_data(point, normal, domain_index, result):
            x, y, z = point - source
            r = np.sqrt(x**2 + y**2 + z**2)
            kr = k * r
            h1kr = (-1j - kr) * np.exp(1j * kr) / (kr * kr)
            scale = h1kr / r
            field = [-y * scale, x * scale, 0.]
            result[:] = np.cross(field, normal)

        def eval_exact_neumann_data(point, normal, domain_index, result):
            x, y, z = point - source
            r = np.sqrt(x**2 + y**2 + z**2)
            kr = k * r
            h1kr = (-1j - kr) * np.exp(1j * kr) / (kr * kr)
            h1kr_deriv = ((1. + 1j - 1j * kr) * (1. + 1j + kr) *
                          np.exp(1j * kr) / (kr * kr * r))
            xy_factor = (h1kr - r * h1kr_deriv) / (r * r * r)
            curl = [x * z * xy_factor,
                    y * z * xy_factor,
                    ((x * x + y * y + 2 * z * z) * h1kr + r *
                     (x * x + y * y) * h1kr_deriv) /
                    (r * r * r)]
            result[:] = np.cross(curl, normal) / (1j * k)

        def eval_exact_solution(point):
            x, y, z = point - source
            r = np.sqrt(x**2 + y**2 + z**2)
            kr = k * r
            h1kr = (-1j - kr) * np.exp(1j * kr) / (kr * kr)
            scale = h1kr / r
            return np.array([-y * scale, x * scale, 0.])

        space = bempp.api.function_space(grid, "RWG", 0)

        dirichlet_data = bempp.api.GridFunction(space, fun=eval_dirichlet_data)
        neumann_data = bempp.api.GridFunction(
            space, fun=eval_exact_neumann_data)

        #pylint: disable=no-member
        eval_point = np.array([[3, 2, 1]]).T

        efie_pot = bempp.api.operators.potential.maxwell.electric_field(
            space, eval_point, k, parameters=parameters)
        mfie_pot = bempp.api.operators.potential.maxwell.magnetic_field(
            space, eval_point, k, parameters=parameters)

        expected = eval_exact_solution(eval_point[:, 0])
        #pylint: disable=unsubscriptable-object
        actual = (-efie_pot * neumann_data - mfie_pot * dirichlet_data)[:, 0]
        rel_error = np.linalg.norm(expected - actual) / np.linalg.norm(actual)

        self.assertTrue(
            rel_error < 1E-3,
            msg="Actual error: {0}".format(rel_error) +
            ". Expected error: 1E-3")

    def test_electric_far_field(self):
        """Test the electric far field operator."""

        # This example computes the far-field pattern for the solution of
        # $$
        # \text{curl}~\text{curl} E - k^2 E = 0
        # $$
        # with boundary data given from an analytic Maxwell solution
        # in the exterior of the sphere as
        # $$
        # E\times n = h_1^{(1)}(kr)e_{\phi}\times n
        # $$
        # Here, $e_{\phi}$ is the unit vector along the $\phi$
        # derivative $\frac{dx}{d\phi}$ in spherical $(r,\theta,\phi)$
        # coordinates.
        #
        # The far-field pattern is given analytically by $-e_{\phi}$.

        parameters = bempp.api.common.global_parameters()
        parameters.assembly.boundary_operator_assembly_type = 'dense'
        parameters.assembly.potential_operator_assembly_type = 'dense'
        grid = bempp.api.shapes.regular_sphere(4)

        rwg_space = bempp.api.function_space(grid, 'RWG', 0)
        snc_space = bempp.api.function_space(grid, "SNC", 0)

        k = 1

        def dirichlet_data(x, n, domain_index, res):
            r = np.sqrt(x[0]**2 + x[1]**2 + x[2]**2)
            kr = k * r
            h1kr = (-1j - kr) * np.exp(1j * kr) / (kr * kr)
            field = h1kr * np.array([-x[1] / r, x[0] / r, 0.])
            res[:] = np.cross(field, n)

        grid_fun = bempp.api.GridFunction(rwg_space, fun=dirichlet_data)

        ident = bempp.api.operators.boundary.sparse.identity(
            rwg_space, rwg_space, snc_space)
        efie = bempp.api.operators.boundary.maxwell.electric_field(
            rwg_space, rwg_space, snc_space, k, parameters=parameters,
            use_projection_spaces=False)

        sol = bempp.api.linalg.lu(efie, ident * grid_fun)

        from bempp.api.operators.far_field.maxwell import \
            electric_field as electric_far_field
        npoints = 100
        theta = np.linspace(0, 2 * np.pi, npoints)
        points = np.vstack([np.cos(theta), np.sin(
            theta), np.zeros(100, dtype='float64')])

        far_field = electric_far_field(
            rwg_space, points, k, parameters=parameters) * sol
        exact = np.vstack([points[1], -points[0], np.zeros(100)])

        rel_error = np.linalg.norm(far_field - exact) / np.linalg.norm(exact)
        self.assertTrue(
            rel_error < 1E-5,
            msg="Actual error: {0}. Expected error: 1E-5".format(rel_error))

if __name__ == "__main__":

    from unittest import main
    main()
