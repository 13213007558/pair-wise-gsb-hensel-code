# frozen_string_literal: true

require "test_helper"

class TestMultiHenselLift < Minitest::Test
  def setup
    @prime = 7
    @g1 = HenselCode::FormalPolynomial.new [1, 1]    # x + 1
    @g2 = HenselCode::FormalPolynomial.new [3, 1]    # x + 3
    @g3 = HenselCode::FormalPolynomial.new [0, 2, 1] # x^2 + 2x
    @factors = [@g1, @g2, @g3]
    # f agrees with the product of the factors modulo 7 but differs over Z,
    # so every lifting level requires a genuine correction
    @f = (@g1 * @g2 * @g3) + (HenselCode::FormalPolynomial.new([2, 0, 1]) * @prime)
  end

  def product(polynomials)
    polynomials.reduce(HenselCode::FormalPolynomial.new([1])) { |acc, g| acc * g }
  end

  def test_initialize
    lift = HenselCode::MultiHenselLift.new @prime, @f, @factors

    assert_equal HenselCode::MultiHenselLift, lift.class
    assert_equal @prime, lift.prime
    assert_equal @f, lift.polynomial
    assert_equal @factors, lift.factors
  end

  def test_synchronous_lift_of_three_factors
    exponent = 4
    result = HenselCode::MultiHenselLift.new(@prime, @f, @factors).lift(exponent)
    modulus = @prime**exponent

    assert_equal 3, result.factors.size
    assert result.factors.all?(&:monic?)
    # the product of the synchronously lifted factors matches f modulo 7^4
    assert_equal @f.mod(modulus), product(result.factors).mod(modulus)
  end

  def test_factor_correspondence_is_preserved
    result = HenselCode::MultiHenselLift.new(@prime, @f, @factors).lift(3)

    assert_equal @factors.map { |g| g.mod(@prime) }, result.input_factors
    result.factors.each_with_index do |lifted, i|
      assert_equal @factors[i].mod(@prime), lifted.mod(@prime)
    end
  end

  def test_per_level_product_residuals
    exponent = 4
    result = HenselCode::MultiHenselLift.new(@prime, @f, @factors).lift(exponent)

    assert_equal (1..exponent).to_a, result.residuals.keys
    result.residuals.each do |k, residual|
      divisible = residual.coefficients.all? { |c| (c % (@prime**k)).zero? }
      assert divisible, "residual at level #{k} not divisible by #{@prime}^#{k}"
    end
    # the residual at the final level is exactly f minus the product of the lifted factors
    assert_equal @f - product(result.factors), result.residuals[exponent]
  end

  def test_lift_to_n_agrees_with_level_by_level_lifting
    exponent = 5
    lift = HenselCode::MultiHenselLift.new @prime, @f, @factors
    full = lift.lift(exponent)

    (1..exponent).each do |k|
      stepwise = lift.lift(k)
      assert_equal full.factors.map { |g| g.mod(@prime**k) }, stepwise.factors
    end
  end

  def test_non_coprime_factors_are_rejected
    f = @g1 * @g1 * @g2

    assert_raises(HenselCode::NonCoprimePolynomialFactors) do
      HenselCode::MultiHenselLift.new @prime, f, [@g1, @g1, @g2]
    end
  end

  def test_factors_not_splitting_the_polynomial_are_rejected
    other = HenselCode::FormalPolynomial.new [4, 1]

    assert_raises(HenselCode::WrongHenselCodeInputType) do
      HenselCode::MultiHenselLift.new @prime, @f, [@g1, @g2, other]
    end
  end

  def test_non_monic_inputs_are_rejected
    not_monic = HenselCode::FormalPolynomial.new [1, 2]
    f = not_monic * @g2

    assert_raises(ArgumentError) { HenselCode::MultiHenselLift.new @prime, f, [not_monic, @g2] }
    assert_raises(ArgumentError) { HenselCode::MultiHenselLift.new @prime, not_monic, [@g1, @g2] }
  end

  def test_invalid_prime_and_exponent_are_rejected
    assert_raises(ArgumentError) { HenselCode::MultiHenselLift.new 1, @f, @factors }
    lift = HenselCode::MultiHenselLift.new @prime, @f, @factors

    assert_raises(ArgumentError) { lift.lift(0) }
    assert_raises(ArgumentError) { lift.lift(1.5) }
  end
end
