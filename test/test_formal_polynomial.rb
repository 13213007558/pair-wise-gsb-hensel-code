# frozen_string_literal: true

require "test_helper"

class TestFormalPolynomial < Minitest::Test
  def test_builds_coefficients_from_padic_digits
    f = HenselCode::FormalPolynomial.from_padic_coefficients(
      5, 4, [[1, 1, 0, 0], [0, 0, 0, 1]]
    )

    assert_equal [6, 125], f.coefficients
    assert_equal 1, f.degree
    assert_equal 1, f.degree_bound
  end

  def test_builds_coefficients_from_hensel_codes
    digits = HenselCode::FinitePadicExpansion.new(5, 4, Rational(6, 1)).hensel_code
    code = HenselCode::TruncatedFinitePadicExpansion.new(5, 4, Rational(3, 2))
    f = HenselCode::FormalPolynomial.from_padic_coefficients(5, 4, [digits, code])

    assert_equal [6, 314], f.coefficients
  end

  def test_multiplication_does_not_truncate_formal_degree
    f = HenselCode::FormalPolynomial.new(5, 2, [1, 2])
    g = HenselCode::FormalPolynomial.new(5, 2, [3, 4])

    assert_equal [3, 10, 8], (f * g).coefficients
  end
end
