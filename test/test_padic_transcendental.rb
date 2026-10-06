# frozen_string_literal: true

require "test_helper"

class TestPAdicTranscendental < Minitest::Test
  EXPANSION_CLASSES = [
    HenselCode::FinitePadicExpansion,
    HenselCode::TruncatedFinitePadicExpansion
  ].freeze

  def test_log_exp_round_trip_for_prime_two
    EXPANSION_CLASSES.each do |klass|
      number = Rational(13, 5)
      x = klass.new(2, 8, number)

      logarithm = x.log(6)
      result = logarithm.exp(6)

      assert_equal 8, logarithm.padic_residue_mod(64)
      assert_equal x.padic_residue_mod(64), result.padic_residue_mod(64)
      assert_equal 6, logarithm.transcendental_metadata[:effective_precision]
      assert_equal 6, result.transcendental_metadata[:effective_precision]
      assert logarithm.transcendental_metadata[:stop_reason].include?(", >= 6")
      assert result.transcendental_metadata[:stop_reason].include?(", >= 6")
    end
  end

  def test_log_exp_round_trip_for_odd_prime
    EXPANSION_CLASSES.each do |klass|
      number = Rational(4, 7)
      x = klass.new(3, 7, number)

      logarithm = x.log(5)
      result = logarithm.exp(5)

      assert_equal 69, logarithm.padic_residue_mod(243)
      assert_equal x.padic_residue_mod(243), result.padic_residue_mod(243)
      assert_equal 1, logarithm.transcendental_metadata[:input_valuation]
      assert logarithm.transcendental_metadata[:denominator_precision_loss].any?(&:positive?)
      assert result.transcendental_metadata[:denominator_precision_loss].any?(&:positive?)
    end
  end

  def test_log_and_exp_fixed_points
    EXPANSION_CLASSES.each do |klass|
      one = klass.new(5, 6, Rational(1, 1))
      zero = klass.new(5, 6, Rational(0, 1))
      logarithm = one.log(6)
      exponential = zero.exp(6)

      assert_equal 0, logarithm.padic_residue_mod(15_625)
      assert_equal 1, exponential.padic_residue_mod(15_625)
      assert logarithm.transcendental_metadata[:stop_reason].include?("congruent to 1")
      assert exponential.transcendental_metadata[:stop_reason].include?("congruent to 0")
    end
  end

  def test_prime_two_rejects_numbers_that_only_look_close_enough
    EXPANSION_CLASSES.each do |klass|
      assert_raises HenselCode::PAdicConvergenceDomainError do
        klass.new(2, 8, Rational(3, 1)).log(6)
      end

      assert_raises HenselCode::PAdicConvergenceDomainError do
        klass.new(2, 8, Rational(5, 3)).log(6)
      end

      assert_raises HenselCode::PAdicConvergenceDomainError do
        klass.new(2, 8, Rational(2, 3)).exp(6)
      end

      assert_raises HenselCode::PAdicConvergenceDomainError do
        klass.new(2, 8, Rational(6, 5)).exp(6)
      end
    end
  end

  def test_odd_prime_domain_rejects_noncanonical_units_and_nonmultiples
    EXPANSION_CLASSES.each do |klass|
      assert_raises HenselCode::PAdicConvergenceDomainError do
        klass.new(3, 7, Rational(5, 4)).log(5)
      end

      assert_raises HenselCode::PAdicConvergenceDomainError do
        klass.new(3, 7, Rational(2, 1)).log(5)
      end

      assert_raises HenselCode::PAdicConvergenceDomainError do
        klass.new(3, 7, Rational(2, 5)).exp(5)
      end
    end
  end

  def test_prime_two_requires_two_digits_to_detect_domain
    EXPANSION_CLASSES.each do |klass|
      assert_raises HenselCode::InsufficientPAdicPrecisionError do
        klass.new(2, 1, Rational(5, 3)).log(1)
      end

      assert_raises HenselCode::InsufficientPAdicPrecisionError do
        klass.new(2, 1, Rational(4, 3)).exp(1)
      end

      assert_raises HenselCode::InsufficientPAdicPrecisionError do
        klass.new(2, 6, Rational(5, 3)).log(7)
      end
    end
  end

  def test_results_are_padic_objects_without_real_float_values
    x = HenselCode::FinitePadicExpansion.new(2, 8, Rational(13, 5))
    logarithm = x.log(6)

    assert_kind_of HenselCode::FinitePadicExpansion, logarithm
    assert_equal 6, logarithm.exponent
    assert logarithm.to_r.numerator.is_a?(Integer)
    assert logarithm.to_r.denominator.is_a?(Integer)
    assert_kind_of Integer, logarithm.transcendental_metadata[:tail_valuation_lower_bound]
    assert_kind_of Integer, logarithm.transcendental_metadata[:first_omitted_index]
  end

  def test_series_implementation_does_not_use_real_log_or_exp
    source = File.read(File.expand_path("../lib/hensel_code/padic_transcendental.rb", __dir__))

    refute_match(/Math\.log/, source)
    refute_match(/Math\.exp/, source)
  end
end
