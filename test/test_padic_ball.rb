# frozen_string_literal: true

require "test_helper"

class TestPAdicBall < Minitest::Test
  def test_exact_and_uncertain_zero_are_distinct
    exact_zero = HenselCode::PAdicBall.new(5, 0, Float::INFINITY, Float::INFINITY)
    uncertain_zero = HenselCode::PAdicBall.new(5, 0, 2, 2)

    assert exact_zero.exact_zero?
    assert exact_zero.contains_zero?
    refute exact_zero.uncertain_zero?
    assert_equal Float::INFINITY, exact_zero.valuation
    assert_equal Float::INFINITY, exact_zero.absolute_precision

    refute uncertain_zero.exact_zero?
    assert uncertain_zero.uncertain_zero?
    assert_equal 2, uncertain_zero.valuation
    assert_equal 2, uncertain_zero.absolute_precision
    assert_equal [0, 25, 50, 75, 100], uncertain_zero.possible_residues(3)
  end

  def test_ball_containment
    broad = HenselCode::PAdicBall.new(5, 0, 1, 1)
    narrow = HenselCode::PAdicBall.new(5, 5, 1, 2)
    exact = HenselCode::PAdicBall.new(5, 30, 1, Float::INFINITY)

    assert broad.cover?(narrow)
    refute narrow.cover?(broad)
    assert narrow.cover?(exact)
    assert broad.cover?(0)
    refute broad.cover?(1)
  end

  def test_canonical_center_and_known_valuation
    divisible = HenselCode::PAdicBall.new(5, 25, 2, 2)
    unit = HenselCode::PAdicBall.new(5, Rational(7, 2), 0, 2)

    assert_equal 0, divisible.center
    assert_equal 2, divisible.valuation
    assert_equal 7 * 13 % 25, unit.center
    assert_equal 0, unit.valuation
  end

  def test_rejects_unprovable_high_valuation
    error = assert_raises(HenselCode::WrongHenselCodeInputType) do
      HenselCode::PAdicBall.new(5, 0, 3, 2)
    end

    assert_match(/valuation 3/, error.message)
  end

  def test_addition_and_subtraction_precision
    precise = HenselCode::PAdicBall.new(5, 1, 0, 3)
    coarser = HenselCode::PAdicBall.new(5, 2, 0, 1)

    assert_equal 1, (precise + coarser).absolute_precision
    assert_equal 1, (precise - coarser).absolute_precision
  end

  def test_identical_ball_difference_is_uncertain_zero
    ball = HenselCode::PAdicBall.new(5, 1, 0, 2)
    difference = ball - ball
    product = difference * difference

    assert_equal 0, difference.center
    assert_equal 2, difference.valuation
    assert difference.contains_zero?
    refute difference.exact_zero?
    assert_equal [0], difference.possible_residues(2)
    refute_equal [0], difference.possible_residues

    assert_equal 0, product.center
    assert_equal 4, product.valuation
    assert_equal 4, product.absolute_precision
    refute product.exact_zero?
  end

  def test_multiplication_precision
    unit = HenselCode::PAdicBall.new(5, 1, 0, 3)
    divisible = HenselCode::PAdicBall.new(5, 10, 1, 4)
    product = unit * divisible

    assert_equal 1, product.valuation
    assert_equal 4, product.absolute_precision
  end

  def test_exact_zero_times_uncertain_ball_is_exact_zero
    exact_zero = HenselCode::PAdicBall.new(5, 0, Float::INFINITY, Float::INFINITY)
    uncertain_zero = HenselCode::PAdicBall.new(5, 0, 2, 2)

    assert (exact_zero * uncertain_zero).exact_zero?
  end

  def test_unit_inverse_precision
    ball = HenselCode::PAdicBall.new(5, 2, 0, 2)
    inverse = ball.inverse

    assert_equal 0, inverse.valuation
    assert_equal 2, inverse.absolute_precision
    assert_equal 13, inverse.center
    assert inverse.cover?(Rational(1, 2))
    assert (ball * inverse).cover?(1)
  end

  def test_nonunit_inverse_has_negative_absolute_precision
    ball = HenselCode::PAdicBall.new(3, 3, 1, 2)
    inverse = ball.inverse

    assert_equal(-1, inverse.valuation)
    assert_equal 0, inverse.absolute_precision
    assert_equal [1], inverse.unit_residues(1)
    assert inverse.cover?(Rational(1, 3))
  end

  def test_inverse_rejects_balls_containing_zero
    [
      HenselCode::PAdicBall.new(5, 0, 2, 2),
      HenselCode::PAdicBall.new(5, 0, Float::INFINITY, Float::INFINITY)
    ].each do |ball|
      error = assert_raises(ZeroDivisionError) { ball.inverse }
      assert_match(/contains zero/, error.message)
    end
  end

  def test_small_residue_balls_cover_every_addition_subtraction_and_multiplication
    [2, 3, 5].each do |prime|
      precision = 2
      balls = (0...(prime**precision)).map do |center|
        HenselCode::PAdicBall.new(prime, center, precision)
      end

      balls.each do |left_ball|
        balls.each do |right_ball|
          %i[+ - *].each do |operation|
            result = left_ball.public_send(operation, right_ball)
            result_exponent = result.absolute_precision + 1
            modulus = prime**result_exponent
            residues = result.possible_residues(result_exponent)

            left_exponent = operation == :* ? result_exponent - right_ball.valuation : result_exponent
            right_exponent = operation == :* ? result_exponent - left_ball.valuation : result_exponent

            left_ball.possible_residues(left_exponent).each do |left|
              right_ball.possible_residues(right_exponent).each do |right|
                value = left.public_send(operation, right) % modulus
                assert_includes residues, value
                assert result.cover?(value % prime**result.absolute_precision)
              end
            end
          end
        end
      end
    end
  end

  def test_small_residue_balls_cover_every_division
    [2, 3, 5].each do |prime|
      precision = 2
      numerators = (0...(prime**precision)).map do |center|
        HenselCode::PAdicBall.new(prime, center, precision)
      end
      denominators = numerators.select { |ball| ball.valuation.zero? }

      numerators.each do |numerator|
        denominators.each do |denominator|
          quotient = numerator / denominator
          result_exponent = quotient.absolute_precision + 1
          modulus = prime**result_exponent
          residues = quotient.possible_residues(result_exponent)

          numerator_exponent = result_exponent - denominator.valuation
          denominator_exponent = result_exponent - numerator.valuation

          numerator.possible_residues(numerator_exponent).each do |left|
            denominator.possible_residues(denominator_exponent).each do |right|
              inverse = denominator.send(:mod_inverse, right, modulus)
              value = (left * inverse) % modulus
              assert_includes residues, value
              assert quotient.cover?(value % prime**quotient.absolute_precision)
            end
          end
        end
      end
    end
  end

  def test_from_truncated
    code = HenselCode::TruncatedFinitePadicExpansion.new(5, 2, 7)
    ball = HenselCode::PAdicBall.from_truncated(code)

    assert_equal 7, ball.center
    assert_equal 0, ball.valuation
    assert_equal 2, ball.absolute_precision

    zero_code = HenselCode::TruncatedFinitePadicExpansion.new(5, 2, 0)
    zero_ball = HenselCode::PAdicBall.from_truncated(zero_code)

    assert zero_ball.contains_zero?
    refute zero_ball.exact_zero?
  end

  def test_operands_must_have_same_prime
    left = HenselCode::PAdicBall.new(3, 1, 0, 2)
    right = HenselCode::PAdicBall.new(5, 1, 0, 2)

    assert_raises(HenselCode::HenselCodesWithDifferentPrimes) { left + right }
  end
end
