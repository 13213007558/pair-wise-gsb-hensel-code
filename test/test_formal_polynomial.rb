# frozen_string_literal: true

require "test_helper"

class TestFormalPolynomial < Minitest::Test
  def setup
    @prime = 257
    @exponent = 3
    @modulus = @prime**@exponent
  end

  def padic(integer, exponent = @exponent)
    HenselCode::TruncatedFinitePadicExpansion.new(@prime, exponent, Rational(integer % (@prime**exponent), 1))
  end

  def test_initialize_with_integer_coefficients
    f = HenselCode::FormalPolynomial.new(@prime, [1, 2, 3])

    assert_equal HenselCode::FormalPolynomial, f.class
    assert_equal @prime, f.prime
    assert_equal [1, 2, 3], f.coefficients
    assert_equal 2, f.degree
  end

  def test_initialize_with_padic_coefficients
    coefficients = [padic(1), padic(2), 3]
    f = HenselCode::FormalPolynomial.new(@prime, coefficients)

    assert_equal coefficients, f.coefficients
  end

  def test_initialize_with_finite_padic_expansion_coefficients
    digits = [3, 1, 2]
    coefficient = HenselCode::FinitePadicExpansion.new(@prime, @exponent, digits)
    f = HenselCode::FormalPolynomial.new(@prime, [coefficient, 1])

    assert_equal 3 + @prime + (2 * (@prime**2)), f.to_s.split(" + ").first.to_i
  end

  def test_precision_is_infinite_for_exact_coefficients
    f = HenselCode::FormalPolynomial.new(@prime, [1, 2, 3])

    assert_equal Float::INFINITY, f.precision
    assert f.exact?
  end

  def test_precision_is_the_smallest_exponent
    f = HenselCode::FormalPolynomial.new(@prime, [padic(1, 4), padic(2, 3), 5])

    assert_equal 3, f.precision
    refute f.exact?
  end

  def test_validations
    assert_raises(ArgumentError) { HenselCode::FormalPolynomial.new(nil, [1]) }
    assert_raises(ArgumentError) { HenselCode::FormalPolynomial.new(@prime, nil) }
    assert_raises(ArgumentError) { HenselCode::FormalPolynomial.new(@prime, "1") }
    assert_raises(ArgumentError) { HenselCode::FormalPolynomial.new(@prime, []) }
    assert_raises(ArgumentError) { HenselCode::FormalPolynomial.new(@prime, [1.5]) }
  end

  def test_coefficient_prime_mismatch
    coefficient = HenselCode::TruncatedFinitePadicExpansion.new(251, @exponent, Rational(1, 1))

    assert_raises(HenselCode::HenselCodesWithDifferentPrimes) do
      HenselCode::FormalPolynomial.new(@prime, [coefficient])
    end
  end

  def test_operand_validations
    f = HenselCode::FormalPolynomial.new(@prime, [1, 2])
    g = HenselCode::FormalPolynomial.new(251, [1, 2])

    assert_raises(HenselCode::WrongHenselCodeInputType) { f.resultant([1, 2]) }
    assert_raises(HenselCode::HenselCodesWithDifferentPrimes) { f.resultant(g) }
  end

  def test_zero_polynomial_is_rejected
    f = HenselCode::FormalPolynomial.new(@prime, [0, 0])
    g = HenselCode::FormalPolynomial.new(@prime, [1, 1])

    assert_raises(HenselCode::Error) { f.resultant(g) }
  end

  def test_to_s
    f = HenselCode::FormalPolynomial.new(@prime, [2, 3, 4])

    assert_equal "2 + 3x + 4x^2", f.to_s
  end

  def test_inspect
    f = HenselCode::FormalPolynomial.new(@prime, [2, 3, 4])

    assert_equal "<FormalPolynomial: 2 + 3x + 4x^2>", f.inspect
  end
end
