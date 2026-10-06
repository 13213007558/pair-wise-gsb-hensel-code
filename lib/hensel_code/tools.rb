# frozen_string_literal: true

module HenselCode
  # Required basics of number generation and modular
  module Tools
    extend self

    def random_integer(bits)
      OpenSSL::BN.rand(bits).to_i
    end

    def random_rational(bits)
      numbers = [random_integer(bits)]
      while numbers.size < 2
        num = random_integer(bits)
        numbers << num if numbers.last.gcd(num) == 1
      end
      Rational(*numbers)
    end

    def random_prime(bits)
      # The OpenSSL library only generates random primes from 16 bits on
      # The 3513 is the first 16-bit number
      # which means there are 3512 primes under 16 bits in length
      # Therefore for allowing values for the parameter 'bits' between 2 and 15,
      # we generate all primes from 2 to 15 bits using Ruby's Prime library
      # and select only those with bit length equal to the value of 'bits'
      # so we sample one element
      # We only do this if the value of 'bits' is less than 16, otherwise we
      # use OpenSSL to generate the prime
      if bits >= 2 && bits < 16
        primes = Prime.first(3512)
        primes.select { |p| p.bit_length == bits }.sample
      elsif bits >= 16
        OpenSSL::BN.generate_prime(bits).to_i
      else
        raise BadBitRangeForRandomPrime, "The bit length must be greater than or equal to 2"
      end
    end

    def random_distinct_numbers(type, quantity, bits)
      numbers = [send("random_#{type}", bits)]
      while numbers.size < quantity
        number = send("random_#{type}", bits)
        numbers << number unless numbers.include?(number)
      end
      numbers
    end

    def eea_core(num1, num2, bound = 0)
      setup = bound.zero? ? [1, 0] : [0, 1]
      x0, x1, y0, y1 = [num1, num2] + setup
      i = 1
      while x1 > bound
        q = x0 / x1
        x0, x1 = x1, x0 - (q * x1)
        y0, y1 = y1, y0 - (q * y1)
        i += 1
      end
      [x0, y0, i, x1, y1]
    end

    def extended_gcd(num1, num2)
      if num1.negative?
        x, y = extended_gcd(-num1, num2)
        [x, -y]
      elsif num2.negative?
        x, y = extended_gcd(num1, -num2)
        [x, y]
      else
        x, y = eea_core(num1, num2)
        [x, y]
      end
    end

    def mod_inverse(num, mod)
      x, y, = extended_gcd(num, mod)
      raise ZeroDivisionError, "#{num} has no inverse modulo #{mod}" unless x == 1

      y % mod
    end

    def crt(moduli, remainders)
      g = moduli.inject(:*)
      result = 0
      moduli.zip(remainders) do |modulus, remainder|
        g_prime = g / modulus
        g_prime_inverse = mod_inverse(g_prime, modulus)
        result += ((g_prime * g_prime_inverse * remainder))
      end
      result % g
    end

    def fault_tolerant_gadic_decode(
      moduli,
      remainders,
      numerator_bound,
      denominator_bound = numerator_bound,
      max_errors: 2
    )
      validate_fault_tolerant_arguments(moduli, remainders, numerator_bound, denominator_bound, max_errors)

      minimum_set_size = moduli.size - max_errors
      candidates = {}

      minimum_set_size.upto(moduli.size) do |set_size|
        moduli.each_index.to_a.combination(set_size) do |indices|
          subset_moduli = indices.map { |index| moduli[index] }
          subset_remainders = indices.map { |index| remainders[index] }
          modulus_product = subset_moduli.inject(:*)
          next unless modulus_product > 2 * numerator_bound * denominator_bound

          rational = bounded_rational_reconstruction(
            subset_moduli,
            subset_remainders,
            numerator_bound,
            denominator_bound
          )
          add_fault_tolerant_candidate(candidates, rational, modulus_product) if rational
        end
      end

      candidates.map do |rational, candidate|
        matching_indices = moduli.each_index.select do |index|
          remainder_matches_rational?(remainders[index], rational, moduli[index])
        end
        matching_product = matching_indices.map { |index| moduli[index] }.inject(:*)

        candidate[:consistent_indices] = matching_indices
        candidate[:inconsistent_indices] = moduli.each_index.to_a - matching_indices
        candidate[:consistent_moduli] = matching_indices.map { |index| moduli[index] }
        candidate[:inconsistent_moduli] = candidate[:inconsistent_indices].map { |index| moduli[index] }
        candidate[:consistent_modulus_product] = matching_product
        candidate[:consistent_modulus_product_condition] = product_condition(
          matching_product,
          2 * numerator_bound * denominator_bound
        )
        candidate[:certified] = matching_product > 2 * numerator_bound * denominator_bound
      end

      certified_candidates = candidates.select { |_, candidate| candidate[:certified] }
      fault_tolerant_result(certified_candidates, numerator_bound, denominator_bound)
    end

    alias fault_tolerant_decode fault_tolerant_gadic_decode

    private

    def validate_fault_tolerant_arguments(moduli, remainders, numerator_bound, denominator_bound, max_errors)
      unless moduli.is_a?(Array) && remainders.is_a?(Array) && moduli.size == remainders.size && !moduli.empty?
        raise ArgumentError, "moduli and remainders must be non-empty arrays of the same size"
      end

      unless moduli.all? { |modulus| modulus.is_a?(Integer) && modulus > 1 }
        raise ArgumentError, "moduli must contain integers greater than 1"
      end

      unless moduli.combination(2).all? { |modulus, other_modulus| modulus.gcd(other_modulus) == 1 }
        raise ArgumentError, "moduli must be pairwise coprime"
      end

      unless remainders.all? { |remainder| remainder.is_a?(Integer) && remainder >= 0 }
        raise ArgumentError, "remainders must contain non-negative integers"
      end

      moduli.zip(remainders).each do |modulus, remainder|
        next if remainder < modulus

        raise ArgumentError, "remainder #{remainder} must be smaller than modulus #{modulus}"
      end

      unless [numerator_bound, denominator_bound].all? { |bound| bound.is_a?(Integer) && bound.positive? }
        raise ArgumentError, "numerator and denominator bounds must be positive integers"
      end

      unless max_errors.is_a?(Integer) && max_errors >= 0
        raise ArgumentError, "max_errors must be a non-negative integer"
      end

      return unless max_errors >= moduli.size

      raise ArgumentError, "max_errors must be smaller than the number of moduli"
    end

    def bounded_rational_reconstruction(moduli, remainders, numerator_bound, denominator_bound)
      modulus_product = moduli.inject(:*)
      encoded_integer = crt(moduli, remainders)
      return Rational(0, 1) if encoded_integer.zero?

      previous_remainder, remainder = modulus_product, encoded_integer
      previous_coefficient, coefficient = 0, 1

      while remainder.positive?
        quotient = previous_remainder / remainder
        previous_remainder, remainder = remainder, previous_remainder - (quotient * remainder)
        previous_coefficient, coefficient = coefficient, previous_coefficient - (quotient * coefficient)
        next unless remainder <= numerator_bound && coefficient.abs <= denominator_bound

        rational = Rational(remainder, coefficient)
        next unless valid_reconstructed_rational?(
          rational,
          encoded_integer,
          modulus_product,
          numerator_bound,
          denominator_bound
        )

        return rational
      end

      nil
    end

    def valid_reconstructed_rational?(rational, encoded_integer, modulus_product, numerator_bound, denominator_bound)
      rational.numerator.abs <= numerator_bound &&
        rational.denominator <= denominator_bound &&
        rational.denominator.gcd(modulus_product) == 1 &&
        remainder_matches_rational?(encoded_integer, rational, modulus_product)
    end

    def remainder_matches_rational?(remainder, rational, modulus)
      (rational.numerator * mod_inverse(rational.denominator, modulus)) % modulus == remainder % modulus
    end

    def add_fault_tolerant_candidate(candidates, rational, modulus_product)
      return if candidates.key?(rational)

      candidates[rational] = {
        rational: rational,
        witness_modulus_product: modulus_product
      }
    end

    def product_condition(product, threshold)
      "#{product} > #{threshold}"
    end

    def fault_tolerant_result(candidates, numerator_bound, denominator_bound)
      threshold = 2 * numerator_bound * denominator_bound
      base_result = {
        numerator_bound: numerator_bound,
        denominator_bound: denominator_bound,
        required_modulus_product: threshold
      }

      if candidates.size == 1
        rational, candidate = candidates.first
        base_result.merge(
          status: :recovered,
          rational: rational,
          consistent_indices: candidate[:consistent_indices],
          inconsistent_indices: candidate[:inconsistent_indices],
          consistent_moduli: candidate[:consistent_moduli],
          inconsistent_moduli: candidate[:inconsistent_moduli],
          consistent_modulus_product: candidate[:consistent_modulus_product],
          consistent_modulus_product_condition: candidate[:consistent_modulus_product_condition]
        )
      elsif candidates.size > 1
        base_result.merge(status: :ambiguous, candidates: candidates.values)
      else
        base_result.merge(status: :unrecoverable)
      end
    end

    def polynomial_variable(index)
      i = index > 1 ? 2 : index
      ["", "p", "p^#{index}"][i]
    end
  end
end
