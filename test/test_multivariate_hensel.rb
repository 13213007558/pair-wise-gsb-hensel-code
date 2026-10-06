# frozen_string_literal: true

require "test_helper"

class TestMultivariateHensel < Minitest::Test
  def setup
    @prime = 5
    # Mutually coupled system: x0 appears in f2 and x1 appears in f1,
    # so neither variable can be lifted by an independent univariate solve.
    # f1 = x0^2 + x1 - 3
    # f2 = x1^2 + x0
    @f1 = HenselCode::MultivariatePolynomial.new(2, [2, 0] => 1, [0, 1] => 1, [0, 0] => -3)
    @f2 = HenselCode::MultivariatePolynomial.new(2, [0, 2] => 1, [1, 0] => 1)
    @equations = [@f1, @f2]
    @seed = [1, 2] # f1(1, 2) = 0, f2(1, 2) = 5, both 0 mod 5
    @exponent = 4
  end

  def test_polynomial_is_structured_not_string
    assert_equal 0, @f1.evaluate(@seed)
    assert_equal 5, @f2.evaluate(@seed)
    assert_equal 0, @f2.evaluate_mod(@seed, @prime)

    df1_dx0 = @f1.partial(0)
    assert_equal({ [1, 0] => 2 }, df1_dx0.terms)
    assert_equal 2, df1_dx0.evaluate(@seed)
  end

  def test_lift_coupled_system
    result = HenselCode::MultivariateHensel.lift(
      equations: @equations, prime: @prime, seed: @seed, exponent: @exponent
    )

    assert result.applicable?
    assert_equal @prime**@exponent, result.modulus
    assert_equal @exponent, result.precision
    assert_equal 2, result.solution.size

    # the lifted solution reduces to the seed modulo p
    assert_equal @seed, result.solution.map { |value| value % @prime }

    # final substitution back into all equations vanishes modulo p^N
    assert_equal [0, 0], result.final_residuals
    assert result.verified?
    @equations.each do |equation|
      assert_equal 0, equation.evaluate_mod(result.solution, @prime**@exponent)
    end
  end

  def test_each_level_reports_residuals_and_precision
    result = HenselCode::MultivariateHensel.lift(
      equations: @equations, prime: @prime, seed: @seed, exponent: @exponent
    )

    assert_equal (@exponent - 1), result.levels.size
    result.levels.each_with_index do |level, index|
      assert_equal index + 1, level.exponent
      assert_equal @prime**(index + 1), level.modulus
      assert_equal index + 2, level.precision
      # residuals at level k are exact and divisible by p^k
      level.residuals.each do |residual|
        assert_equal 0, residual % level.modulus
      end
      assert_equal 2, level.correction.size
    end
  end

  def test_solution_is_genuinely_coupled
    result = HenselCode::MultivariateHensel.lift(
      equations: @equations, prime: @prime, seed: @seed, exponent: @exponent
    )
    x0, x1 = result.solution

    # both unknowns appear in both equations, so the pair cannot be
    # assembled from independent univariate roots
    assert @f1.terms.keys.any? { |exponents| exponents[1].positive? }
    assert @f2.terms.keys.any? { |exponents| exponents[0].positive? }
    # changing one coordinate alone breaks the system modulo p^N
    refute_equal 0, @f1.evaluate_mod([x0 + @prime, x1], @prime**@exponent)
    refute_equal 0, @f2.evaluate_mod([x0, x1 + @prime], @prime**@exponent)
  end

  def test_singular_seed_is_not_applicable
    # f1 = x0^2, f2 = x1^2 - 1 with seed (0, 1):
    # J = [[2*x0, 0], [0, 2*x1]] = [[0, 0], [0, 2]], singular modulo 5
    g1 = HenselCode::MultivariatePolynomial.new(2, [2, 0] => 1)
    g2 = HenselCode::MultivariatePolynomial.new(2, [0, 2] => 1, [0, 0] => -1)

    result = HenselCode::MultivariateHensel.lift(
      equations: [g1, g2], prime: @prime, seed: [0, 1], exponent: @exponent
    )

    refute result.applicable?
    assert_match(/Jacobian is singular/, result.reason)
    assert_nil result.solution
    refute result.verified?
  end

  def test_non_solution_seed_is_not_applicable
    result = HenselCode::MultivariateHensel.lift(
      equations: @equations, prime: @prime, seed: [0, 0], exponent: @exponent
    )

    refute result.applicable?
    assert_match(/not a solution modulo/, result.reason)
  end

  def test_equations_must_be_structured_polynomials
    error = assert_raises(ArgumentError) do
      HenselCode::MultivariateHensel.lift(
        equations: ["x0^2 + x1 - 3", "x1^2 + x0"],
        prime: @prime, seed: @seed, exponent: @exponent
      )
    end
    assert_match(/MultivariatePolynomial/, error.message)
  end

  def test_system_must_be_square
    assert_raises(ArgumentError) do
      HenselCode::MultivariateHensel.lift(
        equations: [@f1], prime: @prime, seed: @seed, exponent: @exponent
      )
    end
  end

  def test_larger_exponent
    exponent = 8
    result = HenselCode::MultivariateHensel.lift(
      equations: @equations, prime: @prime, seed: @seed, exponent: exponent
    )

    assert result.applicable?
    assert_equal exponent, result.precision
    assert_equal exponent - 1, result.levels.size
    @equations.each do |equation|
      assert_equal 0, equation.evaluate_mod(result.solution, @prime**exponent)
    end
  end
end

