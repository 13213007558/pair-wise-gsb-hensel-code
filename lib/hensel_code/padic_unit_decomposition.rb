# frozen_string_literal: true

module HenselCode
  # Multiplicative decomposition of residues modulo p^N
  module PAdicUnitDecomposition
    class << self
      include Tools

      def decompose(prime, exponent, residue, exact_zero_input = false)
        unless Prime.prime?(prime)
          raise ArgumentError, "prime must be a prime number and it was #{prime}"
        end
        raise ArgumentError, "exponent must be positive" unless exponent.positive?

        modulus = prime**exponent
        normalized_residue = residue % modulus
        status = zero_status(normalized_residue, exact_zero_input)
        result = base_result(prime, exponent, modulus, normalized_residue, status, exact_zero_input)
        return result if normalized_residue.zero?

        valuation, unit = valuation_and_unit(prime, exponent, normalized_residue)
        unit_exponent = exponent - valuation
        unit_modulus = prime**unit_exponent
        root = teichmuller_root(prime, unit_exponent, unit)
        principal_unit = unit * mod_inverse(root, unit_modulus) % unit_modulus
        result.merge(
          valuation: valuation,
          valuation_power: prime**valuation,
          unit_exponent: unit_exponent,
          unit_modulus: unit_modulus,
          unit: unit,
          teichmuller_root: root,
          teichmuller_root_modulus: unit_modulus,
          principal_unit: principal_unit,
          principal_unit_storage_modulus: unit_modulus,
          **root_convention(prime, unit_exponent, unit, root)
        )
      end

      private

      def zero_status(residue, exact_zero_input)
        if exact_zero_input
          :exact_zero
        elsif residue.zero?
          :zero_modulo_prime_power
        else
          :nonzero
        end
      end

      def base_result(prime, exponent, modulus, residue, status, exact_zero_input)
        {
          prime: prime,
          exponent: exponent,
          modulus: modulus,
          residue: residue,
          status: status,
          exact_zero: exact_zero_input,
          zero_modulo_prime_power: residue.zero?
        }
      end

      def valuation_and_unit(prime, exponent, residue)
        unit = residue
        valuation = 0
        while valuation < exponent && (unit % prime).zero?
          unit /= prime
          valuation += 1
        end
        [valuation, unit % (prime**(exponent - valuation))]
      end

      def teichmuller_root(prime, unit_exponent, unit)
        return signed_two_adic_root(unit_exponent, unit) if prime == 2

        root = unit % prime
        current_exponent = 1
        while current_exponent < unit_exponent
          next_exponent = current_exponent + 1
          next_modulus = prime**next_exponent
          root = newton_lift(prime, root % next_modulus, next_modulus)
          current_exponent = next_exponent
        end
        root
      end

      def newton_lift(prime, root, modulus)
        root_power_minus_one = modular_power(root, prime - 1, modulus)
        function_value = (root * root_power_minus_one - root) % modulus
        derivative = (prime * root_power_minus_one - 1) % modulus
        (root - (function_value * mod_inverse(derivative, modulus))) % modulus
      end

      def modular_power(base, exponent, modulus)
        result = 1
        base %= modulus
        while exponent.positive?
          result = (result * base) % modulus if exponent.odd?
          base = (base * base) % modulus
          exponent /= 2
        end
        result
      end

      def signed_two_adic_root(unit_exponent, unit)
        return 1 if unit_exponent == 1

        unit_modulus = 2**unit_exponent
        (unit % 4) == 1 ? 1 : unit_modulus - 1
      end

      def root_convention(prime, unit_exponent, unit, root)
        if prime == 2
          two_adic_convention(unit_exponent, unit)
        else
          {
            root_equation_exponent: prime,
            root_equation_target: root,
            root_order_exponent: prime - 1,
            root_order_target: 1,
            principal_unit_modulus: prime,
            sign: nil,
            teichmuller_convention: :standard,
            principal_unit_convention: :one_plus_padic_integers
          }
        end
      end

      def two_adic_convention(unit_exponent, unit)
        if unit_exponent == 1
          {
            root_equation_exponent: 1,
            root_equation_target: 1,
            root_order_exponent: 1,
            root_order_target: 1,
            principal_unit_modulus: 2,
            sign: :positive,
            teichmuller_convention: :unsigned_low_precision,
            principal_unit_convention: :one_plus_two_adic_integers
          }
        else
          {
            root_equation_exponent: 2,
            root_equation_target: 1,
            root_order_exponent: 2,
            root_order_target: 1,
            principal_unit_modulus: 4,
            sign: (unit % 4) == 1 ? :positive : :negative,
            teichmuller_convention: :signed,
            principal_unit_convention: :one_plus_four_two_adic_integers
          }
        end
      end
    end
  end
end
