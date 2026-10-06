# frozen_string_literal: true

require "test_helper"

class TestPAdicBall < Minitest::Test
  def test_exact_zero_is_distinct_from_zero_modulo_p_power
    finite_zero = HenselCode::PAdicBall.new(5, 0, absolute_precision: 3)
    exact_zero = HenselCode::PAdicBall.exact(5, 0)

    assert_equal 0, finite_zero.center
    assert_equal 3, finite_zero.valuation
    assert_equal 3, finite_zero.absolute_precision
    assert finite_zero.contains_zero?
    refute finite_zero.exact_zero?
    assert finite_zero.uncertain?

    assert_equal 0, exact_zero.center
    assert_equal(Float::INFINITY, exact_zero.valuation)
    assert_equal(Float::INFINITY, exact_zero.absolute_precision)
    assert exact_zero.contains_zero?
    assert exact_zero.exact_zero?
    refute exact_zero.uncertain?
  end

  def test_zero_balls_can_have_different_reliable_precision
    less_precise = HenselCode::PAdicBall.zero(5, 1)
    more_precise = HenselCode::PAdicBall.zero(5, 3)

    refute_equal less_precise, more_precise
    assert_equal 1, less_precise.valuation
    assert_equal 3, more_precise.valuation
    assert less_precise.cover?(5)
    refute less_precise.cover?(1)
    assert more_precise.cover?(125)
    refute more_precise.cover?(25)
  end

  def test_constructor_normalizes_centers_to_canonical_residues
    ball = HenselCode::PAdicBall.new(5, -1, absolute_precision: 2)
    rational_ball = HenselCode::PAdicBall.new(5, Rational(1, 5), absolute_precision: 1)
    hidden_zero = HenselCode::PAdicBall.new(5, 25, absolute_precision: 2)
    noncanonical_remnant = HenselCode::PAdicBall.new(5, 25, absolute_precision: 3)

    assert_equal 24, ball.center
    assert_equal 0, ball.valuation
    assert_equal Rational(1, 5), rational_ball.center
    assert_equal(-1, rational_ball.valuation)
    assert_equal 0, hidden_zero.center
    assert_equal 2, hidden_zero.valuation
    assert_equal 25, noncanonical_remnant.center
    assert_equal 2, noncanonical_remnant.valuation
  end

  def test_constructor_rejects_an_unprovable_valuation
    message = "valuation 3 does not match the provable valuation 2"

    error = assert_raises(HenselCode::WrongHenselCodeInputType) do
      HenselCode::PAdicBall.new(5, 25, valuation: 3, absolute_precision: 3)
    end
    assert_equal message, error.message
  end

  def test_exact_and_cancelled_zero_have_different_arithmetic_precision
    prime = 5
    exact_one = HenselCode::PAdicBall.exact(prime, 1)
    finite_minus_one = HenselCode::PAdicBall.new(prime, -1, absolute_precision: 2)
    exact_zero = HenselCode::PAdicBall.exact(prime, 0)
    finite_unit = HenselCode::PAdicBall.new(prime, 1, absolute_precision: 2)

    cancellation = exact_one + finite_minus_one
    exact_product = exact_zero * finite_unit

    assert_equal 0, cancellation.center
    assert_equal 2, cancellation.valuation
    assert_equal 2, cancellation.absolute_precision
    refute cancellation.exact_zero?
    assert exact_product.exact_zero?
  end

  def test_exact_p_factors_change_precision_without_inventing_zero
    prime = 5
    finite_zero = HenselCode::PAdicBall.zero(prime, 2)
    unit = HenselCode::PAdicBall.new(prime, 1, absolute_precision: 2)
    exact_prime = HenselCode::PAdicBall.exact(prime, prime)
    exact_reciprocal = HenselCode::PAdicBall.exact(prime, Rational(1, prime))
    exact_zero = HenselCode::PAdicBall.exact(prime, 0)

    multiplied = finite_zero * exact_prime
    divided = finite_zero / exact_prime
    lifted_unit = unit * exact_reciprocal

    assert_equal 3, multiplied.valuation
    assert_equal 3, multiplied.absolute_precision
    refute multiplied.exact_zero?
    assert_equal 1, divided.valuation
    assert_equal 1, divided.absolute_precision
    refute divided.exact_zero?
    assert_equal Rational(1, 5), lifted_unit.center
    assert_equal(-1, lifted_unit.valuation)
    assert_equal 1, lifted_unit.absolute_precision
    assert exact_zero.exact_zero?
    assert (exact_zero * exact_prime).exact_zero?
  end

  def test_inverse_is_rejected_when_ball_contains_zero
    finite_zero = HenselCode::PAdicBall.new(5, 0, absolute_precision: 2)
    exact_zero = HenselCode::PAdicBall.exact(5, 0)

    assert_raises(ZeroDivisionError) { finite_zero.inverse }
    assert_raises(ZeroDivisionError) { exact_zero.inverse }
    assert_nil finite_zero.safe_inverse
    assert_nil exact_zero.safe_inverse
    assert_raises(ZeroDivisionError) { finite_one / finite_zero }
  end

  def test_inverse_propagates_reliable_precision
    unit = HenselCode::PAdicBall.new(5, 1, absolute_precision: 2)
    divisible = HenselCode::PAdicBall.new(5, 5, absolute_precision: 3)
    fraction = HenselCode::PAdicBall.new(5, Rational(1, 5), absolute_precision: 1)

    assert_equal 1, unit.inverse.center
    assert_equal 2, unit.inverse.absolute_precision
    assert_equal Rational(1, 5), divisible.inverse.center
    assert_equal(-1, divisible.inverse.valuation)
    assert_equal 1, divisible.inverse.absolute_precision
    assert_equal 5, fraction.inverse.center
    assert_equal 1, fraction.inverse.valuation
    assert_equal 3, fraction.inverse.absolute_precision
  end

  def test_small_modulus_enumeration_covers_every_possible_value
    [
      [2, HenselCode::PAdicBall.exact(2, 0), ball(2, 1, 2)],
      [2, HenselCode::PAdicBall.exact(2, 0), ball(2, Rational(1, 2), 1)],
      [2, ball(2, 1, 2), ball(2, 3, 2)],
      [2, ball(2, 1, 2), ball(2, -1, 2)],
      [2, ball(2, 0, 2), ball(2, 1, 2)],
      [2, ball(2, 2, 2), ball(2, 0, 2)],
      [2, ball(2, Rational(1, 2), 1), ball(2, 1, 1)],
      [2, HenselCode::PAdicBall.exact(2, 1), ball(2, -1, 1)],
      [2, HenselCode::PAdicBall.exact(2, 2), ball(2, 1, 2)],
      [2, HenselCode::PAdicBall.exact(2, Rational(1, 2)), ball(2, 1, 1)],
      [2, ball(2, 0, 2), HenselCode::PAdicBall.exact(2, 2)],
      [2, ball(2, 0, 2), HenselCode::PAdicBall.exact(2, Rational(1, 2))],
      [3, ball(3, 1, 1), ball(3, 2, 1)],
      [3, ball(3, 0, 1), ball(3, 1, 1)],
      [3, ball(3, Rational(1, 3), 1), ball(3, 1, 1)]
    ].each do |prime, first, second|
      enumerate_balls(prime, first, second)
    end
  end

  def test_from_truncated_finite_padic_expansion
    encoded = HenselCode::TruncatedFinitePadicExpansion.new(5, 3, 0)
    converted = HenselCode::PAdicBall.from_truncated(encoded)

    assert_equal 0, converted.center
    assert_equal 3, converted.valuation
    assert_equal 3, converted.absolute_precision
    assert converted.uncertain?
  end

  def test_operands_must_use_the_same_prime
    first = HenselCode::PAdicBall.new(5, 1, absolute_precision: 2)
    second = HenselCode::PAdicBall.new(7, 1, absolute_precision: 2)

    assert_raises(HenselCode::HenselCodesWithDifferentPrimes) { first + second }
    assert_raises(HenselCode::IncompatibleOperandTypes) { first + 1 }
  end

  def test_construction_and_modulus_validations
    assert_raises(HenselCode::WrongHenselCodeInputType) do
      HenselCode::PAdicBall.new(4, 1, absolute_precision: 2)
    end
    assert_raises(HenselCode::WrongHenselCodeInputType) do
      HenselCode::PAdicBall.new(5, Complex(2, 3), absolute_precision: 2)
    end
    assert_raises(HenselCode::WrongHenselCodeInputType) do
      HenselCode::PAdicBall.new(5, 1, 2, precision: 3)
    end
    assert_equal 25, HenselCode::PAdicBall.new(5, 1, absolute_precision: 2).modulus
    assert_raises(HenselCode::WrongHenselCodeInputType) do
      HenselCode::PAdicBall.exact(5, 1).modulus
    end
    assert_raises(HenselCode::WrongHenselCodeInputType) do
      HenselCode::PAdicBall.new(5, Rational(1, 5), absolute_precision: -1).modulus
    end
  end

  private

  def finite_one
    HenselCode::PAdicBall.new(5, 1, absolute_precision: 2)
  end

  def ball(prime, center, precision)
    HenselCode::PAdicBall.new(prime, center, absolute_precision: precision)
  end

  def enumerate_balls(prime, first, second)
    finite_precisions = [first, second].reject(&:exact?).map(&:absolute_precision)
    fine_precision = finite_precisions.max + 2
    shift = -[first.valuation, second.valuation].min
    shift = 0 if shift.negative?
    modulus = prime**(fine_precision + shift)

    first_values = residues(prime, first, shift, fine_precision, modulus)
    second_values = residues(prime, second, shift, fine_precision, modulus)

    refute_empty first_values
    refute_empty second_values

    first_values.each do |first_residue|
      second_values.each do |second_residue|
        first_value = Rational(first_residue, prime**shift)
        second_value = Rational(second_residue, prime**shift)

        %i[+ - *].each do |operation|
          result = first.send(operation, second)
          actual = first_value.send(operation, second_value)
          assert result.cover?(actual), failure_message(result, actual)
        end

        next if second.contains_zero?

        quotient = first / second
        inverse = second.inverse
        actual_quotient = first_value / second_value
        assert quotient.cover?(actual_quotient), failure_message(quotient, actual_quotient)
        assert inverse.cover?(1 / second_value), failure_message(inverse, 1 / second_value)
      end
    end
  end

  def residues(prime, ball, shift, fine_precision, modulus)
    scaled_center = ball.center * prime**shift
    return [scaled_center.to_i % modulus] if ball.exact?

    variation_exponent = shift + ball.absolute_precision
    return [scaled_center.to_i % modulus] if variation_exponent >= fine_precision + shift

    step = prime**variation_exponent
    count = modulus / step
    base = scaled_center % modulus
    Array.new(count) { |index| (base + (index * step)) % modulus }
  end

  def failure_message(ball, actual)
    "#{ball.inspect} does not contain #{actual}"
  end
end
