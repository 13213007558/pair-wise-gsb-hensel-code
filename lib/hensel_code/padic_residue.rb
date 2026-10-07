# frozen_string_literal: true

module HenselCode
  # Residue of an integer modulo p^precision together with its observed p-adic valuation.
  class PadicResidue
    attr_reader :prime, :precision, :residue, :valuation

    def initialize(prime, precision, residue)
      validate_base!(prime, precision)
      @prime = prime
      @precision = precision
      @residue = residue % modulus
      @valuation = compute_valuation
    end

    def modulus
      prime**precision
    end

    def valuation_determined?
      valuation < precision
    end

    def valuation_uncertain?
      !valuation_determined?
    end

    def known_precision
      valuation_determined? ? precision - valuation : 0
    end

    def unresolved_precision
      valuation_determined? ? 0 : precision
    end

    def unit_residue
      return nil unless valuation_determined?

      residue / (prime**valuation)
    end

    def zero?
      residue.zero?
    end

    def to_h
      {
        prime: prime,
        precision: precision,
        residue: residue,
        valuation: valuation_determined? ? valuation : nil,
        valuation_lower_bound: valuation,
        valuation_determined: valuation_determined?,
        known_unit_precision: known_precision,
        unresolved_precision: unresolved_precision
      }
    end

    def to_s
      residue.to_s
    end

    def inspect
      format("<PadicResidue: %d mod %d^%d, valuation=%s>", residue, prime, precision,
             valuation_determined? ? valuation : ">= #{precision}")
    end

    private

    def validate_base!(prime, precision)
      raise ArgumentError, "prime must be an integer" unless prime.is_a?(Integer)
      raise ArgumentError, "prime must be positive" unless prime.positive?
      raise ArgumentError, "precision must be positive" unless precision.positive?
    end

    def compute_valuation
      value = residue
      exponent = 0
      while exponent < precision && (value % prime).zero?
        value /= prime
        exponent += 1
      end
      exponent
    end
  end
end
