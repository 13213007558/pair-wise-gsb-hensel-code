# frozen_string_literal: true

module HenselCode
  # Convergent p-adic power series computed directly in Z/p^nZ.
  module PAdicSeries
    Result = Struct.new(
      :prime,
      :function,
      :target_precision,
      :effective_precision,
      :value,
      :terms_evaluated,
      :first_omitted_term,
      :first_omitted_valuation,
      :stop_reason,
      keyword_init: true
    )

    def log(target_precision)
      prime, precision = series_input!(target_precision)
      input_modulus = prime**exponent
      residue = padic_series_residue

      unless (residue % prime).nonzero?
        raise ConvergenceDomainError, log_domain_message(prime)
      end

      argument = (residue - 1) % input_modulus
      return zero_series_result(:log, prime, precision, 0) if argument.zero?

      valuation, unit = padic_parts(argument, prime)
      ensure_convergence!(:log, prime, valuation)

      first_omitted = cutoff_index(prime, precision, valuation)
      value = log_series(prime, precision, unit, valuation, first_omitted)
      omitted_valuation = (valuation * first_omitted) - integer_valuation(first_omitted, prime)
      series_result(
        function: :log,
        prime: prime,
        precision: precision,
        value: value,
        first_omitted: first_omitted,
        omitted_valuation: omitted_valuation,
        terms_evaluated: first_omitted - 1
      )
    end

    def exp(target_precision)
      prime, precision = series_input!(target_precision)
      argument = padic_series_residue % (prime**exponent)

      return zero_series_result(:exp, prime, precision, 1) if argument.zero?

      valuation, unit = padic_parts(argument, prime)
      ensure_convergence!(:exp, prime, valuation)

      first_omitted = cutoff_index(prime, precision, valuation)
      value = exp_series(prime, precision, unit, valuation, first_omitted)
      omitted_valuation = (valuation * first_omitted) - factorial_valuation(first_omitted, prime)
      series_result(
        function: :exp,
        prime: prime,
        precision: precision,
        value: value,
        first_omitted: first_omitted,
        omitted_valuation: omitted_valuation,
        terms_evaluated: first_omitted
      )
    end

    private

    def series_input!(target_precision)
      unless target_precision.is_a?(Integer) && target_precision.positive?
        raise ArgumentError, "target precision must be a positive integer"
      end
      unless Prime.prime?(prime)
        raise ArgumentError, "p-adic series require a prime, got #{prime}"
      end

      required = [target_precision, convergence_radius(prime)].max
      if exponent < required
        raise InsufficientPrecision,
              "need exponent #{required} to establish #{convergence_radius(prime)} " \
              "initial p-adic digits and target precision #{target_precision}, got #{exponent}"
      end

      [prime, target_precision]
    end

    def padic_series_residue
      if hensel_code.is_a?(Array)
        hensel_code.each_with_index.sum { |digit, index| digit * (prime**index) }
      else
        hensel_code
      end
    end

    def padic_parts(residue, prime)
      valuation = 0
      while (residue % prime).zero?
        residue /= prime
        valuation += 1
      end
      [valuation, residue]
    end

    def convergence_radius(prime)
      prime == 2 ? 2 : 1
    end

    def ensure_convergence!(function, prime, valuation)
      radius = convergence_radius(prime)
      return if valuation >= radius

      raise ConvergenceDomainError, send("#{function}_domain_message", prime)
    end

    def log_domain_message(prime)
      radius = convergence_radius(prime)
      "p-adic log requires x in 1 + #{prime}^#{radius}Z_#{prime}; " \
        "for p=2 use 1+4Z_2, not 1+2Z_2"
    end

    def exp_domain_message(prime)
      radius = convergence_radius(prime)
      "p-adic exp requires y in #{prime}^#{radius}Z_#{prime}; " \
        "for p=2 use 4Z_2, not 2Z_2"
    end

    def cutoff_index(prime, precision, valuation)
      denominator = (valuation * (prime - 1)) - 1
      Rational((precision * (prime - 1)) - 1, denominator).ceil
    end

    def log_series(prime, precision, unit, valuation, first_omitted)
      modulus = prime**precision
      sum = 0
      full_unit_power = 1

      (1...first_omitted).each do |index|
        full_unit_power = (full_unit_power * unit) % modulus
        denominator, denominator_valuation = strip_prime(index, prime)
        shift = (valuation * index) - denominator_valuation
        next if shift >= precision

        term_modulus = prime**(precision - shift)
        term = (full_unit_power % term_modulus) * mod_inverse(denominator, term_modulus)
        term = (-term) % term_modulus if index.even?
        sum = (sum + (term * (prime**shift))) % modulus
      end

      sum
    end

    def exp_series(prime, precision, unit, valuation, first_omitted)
      modulus = prime**precision
      sum = 1
      full_unit_power = 1
      full_factorial_unit = 1
      factorial_valuation = 0

      (1...first_omitted).each do |index|
        full_unit_power = (full_unit_power * unit) % modulus
        factor, factor_valuation = strip_prime(index, prime)
        full_factorial_unit = (full_factorial_unit * factor) % modulus
        factorial_valuation += factor_valuation
        shift = (valuation * index) - factorial_valuation
        next if shift >= precision

        term_modulus = prime**(precision - shift)
        term = (full_unit_power % term_modulus) *
               mod_inverse(full_factorial_unit % term_modulus, term_modulus)
        sum = (sum + (term * (prime**shift))) % modulus
      end

      sum
    end

    def strip_prime(number, prime)
      valuation = 0
      while (number % prime).zero?
        number /= prime
        valuation += 1
      end
      [number, valuation]
    end

    def integer_valuation(number, prime)
      strip_prime(number, prime).last
    end

    def factorial_valuation(number, prime)
      valuation = 0
      while number.positive?
        number /= prime
        valuation += number
      end
      valuation
    end

    def series_result(attributes)
      precision = attributes.fetch(:precision)
      first_omitted = attributes.fetch(:first_omitted)
      omitted_valuation = attributes.fetch(:omitted_valuation)
      stop_reason =
        "tail starts at k=#{first_omitted}; v_p(term)=#{omitted_valuation} >= #{precision}, " \
        "and v_p(k) or v_p(k!) <= (k-1)/(p-1) bounds every later term"

      Result.new(
        prime: attributes.fetch(:prime),
        function: attributes.fetch(:function),
        target_precision: precision,
        effective_precision: precision,
        value: residue_to_padic(precision, attributes.fetch(:value)),
        terms_evaluated: attributes.fetch(:terms_evaluated),
        first_omitted_term: first_omitted,
        first_omitted_valuation: omitted_valuation,
        stop_reason: stop_reason
      )
    end

    def zero_series_result(function, prime, precision, residue)
      constant_count = function == :exp ? 1 : 0
      stop_reason =
        "argument is congruent to 0 modulo #{prime}^#{exponent}; " \
        "all nonconstant series terms already have valuation at least #{exponent}"

      Result.new(
        prime: prime,
        function: function,
        target_precision: precision,
        effective_precision: precision,
        value: residue_to_padic(precision, residue),
        terms_evaluated: constant_count,
        first_omitted_term: 1,
        first_omitted_valuation: exponent,
        stop_reason: stop_reason
      )
    end

    def residue_to_padic(precision, residue)
      if hensel_code.is_a?(Array)
        digits = Array.new(precision) { |index| (residue / (prime**index)) % prime }
        return build_finite_padic(precision, digits)
      end

      build_truncated_padic(precision, residue % (prime**precision))
    end

    def build_finite_padic(precision, digits)
      object = self.class.allocate
      object.instance_variable_set(:@prime, prime)
      object.instance_variable_set(:@exponent, precision)
      object.instance_variable_set(:@n, Integer.sqrt(((prime**precision) - 1) / 2))
      object.instance_variable_set(:@hensel_code, digits)
      object.instance_variable_set(:@polynomial, Polynomial.new(prime, digits))
      object
    end

    def build_truncated_padic(precision, residue)
      object = self.class.allocate
      object.instance_variable_set(:@prime, prime)
      object.instance_variable_set(:@exponent, precision)
      object.instance_variable_set(:@n, Integer.sqrt(((prime**precision) - 1) / 2))
      object.instance_variable_set(:@hensel_code, residue)
      object
    end
  end

  class ConvergenceDomainError < Error; end
  class InsufficientPrecision < Error; end
end
