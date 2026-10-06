# frozen_string_literal: true

module HenselCode
  # p-adic logarithm and exponential power series.
  module PAdicTranscendental
    def log(target_precision = exponent)
      target = validate_transcendental_precision(target_precision, prime == 2 ? 2 : 1)
      modulus = prime**target
      input = (padic_residue_mod(modulus) - 1) % modulus
      validate_log_domain

      if input.zero?
        result = build_padic_result(target, 0)
        metadata = transcendent_metadata(:log, target, input, target, [])
        metadata[:stop_reason] = "input is congruent to 1 modulo #{modulus}"
        metadata[:first_omitted_index] = 1
        metadata[:tail_valuation_lower_bound] = target
        attach_transcendental_metadata(result, metadata)
        return result
      end

      valuation = padic_int_valuation(input, modulus)
      input_unit = input / (prime**valuation)
      terms = []
      sum = 0
      k = 1

      loop do
        denominator_unit, denominator_loss = factor_out_prime(k)
        term_valuation = (k * valuation) - denominator_loss
        term = if term_valuation >= target
                 0
               else
                 unit_modulus = prime**(target - term_valuation)
                 unit_power = input_unit.pow(k, unit_modulus)
                 component = unit_power * mod_inverse(denominator_unit, unit_modulus)
                 component = -component if k.even?
                 (prime**term_valuation) * (component % unit_modulus)
               end
        sum = (sum + term) % modulus
        terms << [term_valuation, denominator_loss]

        tail_bound = (k * valuation) - (k / (prime - 1))
        k += 1
        next unless tail_bound >= target

        result = build_padic_result(target, sum)
        metadata = transcendent_metadata(:log, target, input, valuation, terms)
        metadata[:first_omitted_index] = k
        metadata[:tail_valuation_lower_bound] = tail_bound
        metadata[:stop_reason] = "first omitted term has valuation lower bound #{tail_bound}, >= #{target}"
        attach_transcendental_metadata(result, metadata)
        return result
      end
    end

    def exp(target_precision = exponent)
      target = validate_transcendental_precision(target_precision, prime == 2 ? 2 : 1)
      modulus = prime**target
      input = padic_residue_mod(modulus)
      validate_exp_domain

      if input.zero?
        result = build_padic_result(target, 1)
        metadata = transcendent_metadata(:exp, target, input, target, [[0, 0]])
        metadata[:stop_reason] = "input is congruent to 0 modulo #{modulus}"
        metadata[:first_omitted_index] = 1
        metadata[:tail_valuation_lower_bound] = target
        attach_transcendental_metadata(result, metadata)
        return result
      end

      valuation = padic_int_valuation(input, modulus)
      input_unit = input / (prime**valuation)
      terms = [[0, 0]]
      sum = 1
      factorial_unit = 1
      factorial_p = 0
      k = 1

      loop do
        k_unit, k_p = factor_out_prime(k)
        factorial_unit = (factorial_unit * k_unit) % modulus
        factorial_p += k_p
        term_valuation = (k * valuation) - factorial_p
        term = if term_valuation >= target
                 0
               else
                 unit_modulus = prime**(target - term_valuation)
                 unit_power = input_unit.pow(k, unit_modulus)
                 component = unit_power * mod_inverse(factorial_unit % unit_modulus, unit_modulus)
                 (prime**term_valuation) * (component % unit_modulus)
               end
        sum = (sum + term) % modulus
        terms << [term_valuation, factorial_p]

        tail_bound = (k * valuation) - (k / (prime - 1))
        k += 1
        next unless tail_bound >= target

        result = build_padic_result(target, sum)
        metadata = transcendent_metadata(:exp, target, input, valuation, terms)
        metadata[:first_omitted_index] = k
        metadata[:tail_valuation_lower_bound] = tail_bound
        metadata[:stop_reason] = "first omitted term has valuation lower bound #{tail_bound}, >= #{target}"
        attach_transcendental_metadata(result, metadata)
        return result
      end
    end

    def transcendental_metadata
      @transcendental_metadata
    end

    private

    def validate_transcendental_precision(target, required_domain_precision)
      unless target.is_a?(Integer) && target.positive?
        raise ArgumentError, "target precision must be a positive Integer"
      end

      if target > exponent
        message = "target precision #{target} exceeds input precision #{exponent}"
        raise InsufficientPAdicPrecisionError, message
      end

      if exponent < required_domain_precision
        modulus = prime**required_domain_precision
        message = "at least precision #{required_domain_precision} is required to test convergence modulo #{modulus}"
        raise InsufficientPAdicPrecisionError, message
      end

      target
    end

    def validate_log_domain
      required = prime == 2 ? 2 : 1
      modulus = prime**required
      domain_residue = padic_residue_mod(modulus) % modulus
      return if domain_residue == 1

      domain = prime == 2 ? "x \\equiv 1 (mod 4)" : "x \\equiv 1 (mod #{prime})"
      message = "p-adic log requires #{domain}; got residue #{domain_residue} modulo #{modulus}"
      raise PAdicConvergenceDomainError, message
    end

    def validate_exp_domain
      required = prime == 2 ? 2 : 1
      modulus = prime**required
      domain_residue = padic_residue_mod(modulus) % modulus
      return if domain_residue.zero?

      domain = prime == 2 ? "x \\equiv 0 (mod 4)" : "x \\equiv 0 (mod #{prime})"
      message = "p-adic exp requires #{domain}; got residue #{domain_residue} modulo #{modulus}"
      raise PAdicConvergenceDomainError, message
    end

    def padic_int_valuation(number, modulus)
      valuation = 0
      while valuation < exponent && number % (prime**(valuation + 1)) == 0
        valuation += 1
      end
      [valuation, exponent].min
    end

    def factor_out_prime(number)
      unit = number
      power = 0
      while (unit % prime).zero?
        unit /= prime
        power += 1
      end
      [unit, power]
    end

    def build_padic_result(target, residue)
      self.class.new(prime, target, Rational(residue, 1))
    end

    def transcendent_metadata(function, target, input, valuation, terms)
      {
        function: function,
        prime: prime,
        input_precision: exponent,
        target_precision: target,
        effective_precision: target,
        input_residue: input,
        input_valuation: valuation,
        terms_evaluated: terms.size,
        term_valuations: terms.map(&:first),
        denominator_precision_loss: terms.map(&:last)
      }
    end

    def attach_transcendental_metadata(result, metadata)
      result.instance_variable_set(:@transcendental_metadata, metadata)
      result
    end
  end
end
