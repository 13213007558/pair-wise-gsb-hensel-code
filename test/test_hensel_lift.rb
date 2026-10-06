# frozen_string_literal: true

require "test_helper"

class TestHenselLift < Minitest::Test
  def setup
    @prime = 5
    # f = (x - 1)(x - 2)(x - 3) = x^3 - 6x^2 + 11x - 6
    @polynomial = HenselCode::FormalPolynomial.new([-6, 11, -6, 1])
    @factors = [[-1, 1], [-2, 1], [-3, 1]].map { |c| HenselCode::FormalPolynomial.new(c) }
  end

  def product_of(factors)
    factors.reduce(HenselCode::FormalPolynomial.one) { |acc, factor| acc * factor }
  end

  def test_three_factor_simultaneous_lift
    target_exponent = 4
    lift = HenselCode::HenselLift.lift(@polynomial, @factors, @prime, target_exponent)

    assert_equal target_exponent, lift.levels.size
    lift.levels.each do |level|
      assert_equal @prime**level.exponent, level.modulus
      # per-level residual f - prod(g_i) vanishes modulo the level modulus
      assert(level.residual.coefficients.all? { |c| (c % level.modulus).zero? })
      # lifted factors stay monic
      assert(level.factors.all? { |factor| factor.mod(level.modulus).monic? })
    end

    final_modulus = @prime**target_exponent
    assert (product_of(lift.factors) - @polynomial).mod(final_modulus).zero?

    # the correspondence preserves the order of the input factors and each
    # lifted factor reduces to its own input factor modulo p
    assert_equal @factors.size, lift.correspondence.size
    lift.correspondence.each_with_index do |(input_factor, lifted_factor), i|
      assert_equal @factors[i].mod(@prime), input_factor
      assert_equal input_factor, lifted_factor.mod(@prime)
    end
  end

  def test_lift_to_n_matches_stepwise_lifting
    direct = HenselCode::HenselLift.lift(@polynomial, @factors, @prime, 4)
    halfway = HenselCode::HenselLift.lift(@polynomial, @factors, @prime, 2)
    continued = HenselCode::HenselLift.lift(@polynomial, halfway.factors, @prime, 4, start_exponent: 2)

    assert_equal direct.factors, continued.factors
    continued.levels.each do |level|
      direct_level = direct.levels[level.exponent - 1]

      assert_equal direct_level.factors, level.factors
      assert_equal direct_level.residual, level.residual
    end
  end

  def test_factors_need_not_factor_over_integers
    # f = (x^2 + 1)(x - 1) = x^3 - x^2 + x - 1; x^2 + 1 is irreducible over
    # Z, but x^2 + 1 = (x - 2)(x - 3) (mod 5), so three coprime factors mod 5
    polynomial = HenselCode::FormalPolynomial.new([-1, 1, -1, 1])
    lift = HenselCode::HenselLift.lift(polynomial, @factors, @prime, 3)

    assert (product_of(lift.factors) - polynomial).mod(@prime**3).zero?
    lift.correspondence.each_with_index do |(input_factor, lifted_factor), i|
      assert_equal @factors[i].mod(@prime), input_factor
      assert_equal input_factor, lifted_factor.mod(@prime)
    end
  end

  def test_non_coprime_factors_rejected
    # (x - 1)^2 = x^2 - 2x + 1 with two identical factors mod 5
    polynomial = HenselCode::FormalPolynomial.new([1, -2, 1])
    factors = [[-1, 1], [-1, 1]].map { |c| HenselCode::FormalPolynomial.new(c) }

    assert_raises(HenselCode::NonCoprimeFactors) do
      HenselCode::HenselLift.lift(polynomial, factors, @prime, 3)
    end
  end

  def test_factors_coprime_over_integers_but_not_mod_p_rejected
    # (x - 1)(x - 6) = x^2 - 7x + 6; the factors are distinct over Z but
    # coincide modulo 5
    polynomial = HenselCode::FormalPolynomial.new([6, -7, 1])
    factors = [[-1, 1], [-6, 1]].map { |c| HenselCode::FormalPolynomial.new(c) }

    assert_raises(HenselCode::NonCoprimeFactors) do
      HenselCode::HenselLift.lift(polynomial, factors, @prime, 3)
    end
  end

  def test_product_mismatch_rejected
    factors = [[-1, 1], [-2, 1], [-4, 1]].map { |c| HenselCode::FormalPolynomial.new(c) }

    assert_raises(ArgumentError) do
      HenselCode::HenselLift.lift(@polynomial, factors, @prime, 3)
    end
  end

  def test_non_monic_polynomial_rejected
    polynomial = HenselCode::FormalPolynomial.new([-6, 11, -6, 2])

    assert_raises(ArgumentError) do
      HenselCode::HenselLift.lift(polynomial, @factors, @prime, 3)
    end
  end

  def test_non_prime_modulus_rejected
    assert_raises(ArgumentError) do
      HenselCode::HenselLift.lift(@polynomial, @factors, 4, 3)
    end
  end

  def test_single_factor_lifts_to_the_polynomial_itself
    lift = HenselCode::HenselLift.lift(@polynomial, [@polynomial], @prime, 3)

    assert_equal [@polynomial.mod(@prime**3)], lift.factors
    assert lift.levels.last.residual.mod(@prime**3).zero?
  end
end
