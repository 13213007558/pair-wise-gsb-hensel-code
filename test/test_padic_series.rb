# frozen_string_literal: true

require "test_helper"

class TestPAdicSeries < Minitest::Test
  PADIC_CLASSES = [
    HenselCode::TruncatedFinitePadicExpansion,
    HenselCode::FinitePadicExpansion
  ].freeze

  ROUND_TRIP_CASES = [
    [2, 5, Rational(5, 1)],
    [2, 5, Rational(7, 3)],
    [3, 6, Rational(10, 7)],
    [5, 5, Rational(11, 6)],
    [7, 5, Rational(8, 1)]
  ].freeze

  def test_exp_log_round_trip_for_both_expansions
    PADIC_CLASSES.each do |klass|
      ROUND_TRIP_CASES.each do |prime, precision, rational|
        original = klass.new(prime, precision, rational)
        log_result = original.log(precision)
        exp_result = log_result.value.exp(precision)
        expected = klass.new(prime, precision, rational)

        assert_equal padic_code(expected), padic_code(exp_result.value)
        assert_equal precision, log_result.effective_precision
        assert_equal precision, exp_result.effective_precision
        assert_operator log_result.first_omitted_valuation, :>=, precision
        assert_operator exp_result.first_omitted_valuation, :>=, precision
        assert_match(/tail starts at k=\d+/, log_result.stop_reason)
        assert_match(/v_p\(term\)=\d+ >= #{precision}/, exp_result.stop_reason)
      end
    end
  end

  def test_known_series_values_are_computed_as_padic_residues
    assert_equal 555, HenselCode::TruncatedFinitePadicExpansion.new(5, 4, 6).log(4).value.to_i
    assert_equal 1805, HenselCode::TruncatedFinitePadicExpansion.new(5, 6, 6).log(5).value.to_i
    assert_equal 456, HenselCode::TruncatedFinitePadicExpansion.new(5, 4, 5).exp(4).value.to_i
    assert_equal 2956, HenselCode::TruncatedFinitePadicExpansion.new(5, 6, 5).exp(5).value.to_i
    assert_equal 12, HenselCode::TruncatedFinitePadicExpansion.new(2, 4, 5).log(4).value.to_i
    assert_equal 13, HenselCode::TruncatedFinitePadicExpansion.new(2, 4, Rational(4, 1)).exp(4).value.to_i
  end

  def test_finite_expansion_returns_digits_at_target_precision
    result = HenselCode::FinitePadicExpansion.new(2, 5, Rational(5, 1)).log(4)

    assert_equal [0, 0, 1, 1], result.value.to_a
    assert_equal 4, result.effective_precision
  end

  def test_log_rejects_values_outside_one_plus_p_z_p
    assert_raises HenselCode::ConvergenceDomainError do
      HenselCode::TruncatedFinitePadicExpansion.new(2, 4, 3).log(4)
    end
    assert_raises HenselCode::ConvergenceDomainError do
      HenselCode::TruncatedFinitePadicExpansion.new(2, 4, Rational(5, 3)).log(4)
    end
    assert_raises HenselCode::ConvergenceDomainError do
      HenselCode::TruncatedFinitePadicExpansion.new(3, 4, Rational(5, 4)).log(4)
    end
    assert_raises HenselCode::ConvergenceDomainError do
      HenselCode::TruncatedFinitePadicExpansion.new(5, 4, Rational(7, 6)).log(4)
    end
  end

  def test_exp_rejects_values_outside_p_z_p
    assert_raises HenselCode::ConvergenceDomainError do
      HenselCode::TruncatedFinitePadicExpansion.new(2, 4, 2).exp(4)
    end
    assert_raises HenselCode::ConvergenceDomainError do
      HenselCode::TruncatedFinitePadicExpansion.new(2, 4, Rational(2, 3)).exp(4)
    end
    assert_raises HenselCode::ConvergenceDomainError do
      HenselCode::TruncatedFinitePadicExpansion.new(3, 4, Rational(1, 4)).exp(4)
    end
    assert_raises HenselCode::ConvergenceDomainError do
      HenselCode::TruncatedFinitePadicExpansion.new(5, 4, Rational(1, 6)).exp(4)
    end
  end

  def test_p_two_domain_error_names_four_z_two
    error = assert_raises HenselCode::ConvergenceDomainError do
      HenselCode::TruncatedFinitePadicExpansion.new(2, 4, 3).log(4)
    end

    assert_match(/1 \+ 2\^2Z_2/, error.message)
    assert_match(/1\+4Z_2/, error.message)
  end

  def test_insufficient_expansion_precision_is_rejected
    log_residue = truncated_from_residue(2, 1, 3)
    assert_raises HenselCode::InsufficientPrecision do
      log_residue.log(1)
    end
    exp_residue = truncated_from_residue(2, 1, 2)
    assert_raises HenselCode::InsufficientPrecision do
      exp_residue.exp(1)
    end
    assert_raises HenselCode::InsufficientPrecision do
      HenselCode::TruncatedFinitePadicExpansion.new(5, 3, 6).log(4)
    end
  end

  def test_target_precision_must_be_positive_integer
    h = HenselCode::TruncatedFinitePadicExpansion.new(5, 4, 6)

    assert_raises(ArgumentError) { h.log(0) }
    assert_raises(ArgumentError) { h.exp(-1) }
    assert_raises(ArgumentError) { h.log(4.0) }
  end

  def test_series_are_computed_without_real_transcendentals
    source = File.read(File.expand_path("../lib/hensel_code/padic_series.rb", __dir__))

    refute_match(/\bMath\.(log|exp)\b/, source)
    refute_match(/\bFloat\b/, source)
  end

  private

  def padic_code(padic)
    padic.respond_to?(:to_i) ? padic.to_i : padic.to_a
  end

  def truncated_from_residue(prime, exponent, residue)
    padic = HenselCode::TruncatedFinitePadicExpansion.allocate
    padic.instance_variable_set(:@prime, prime)
    padic.instance_variable_set(:@exponent, exponent)
    padic.instance_variable_set(:@n, Integer.sqrt(((prime**exponent) - 1) / 2))
    padic.instance_variable_set(:@hensel_code, residue)
    padic
  end
end
