# frozen_string_literal: true

module HenselCode
  # Resultant of two formal polynomials together with its p-adic precision
  # analysis.
  #
  # When every coefficient of both polynomials is an exact Integer, the
  # precision is infinite and the value is the exact integer resultant.
  # Otherwise the value is only known modulo p^r, where r is the smallest
  # precision among the p-adic coefficients. In that case a zero residue does
  # NOT mean the resultant is zero: it only means its p-adic valuation is at
  # least r, so no conclusion about common factors can be drawn.
  class ResultantResult
    include HenselCode::Tools

    attr_reader :prime, :precision, :value

    def initialize(prime:, precision:, value:)
      @prime = prime
      @precision = precision
      @value = value
    end

    def exact?
      precision.infinite?
    end

    def modulus
      exact? ? nil : prime**precision
    end

    # canonical representative of the resultant: the exact integer when the
    # precision is infinite, otherwise the residue modulo p^precision
    def to_i
      exact? ? value : value % modulus
    end

    # true when the p-adic valuation of the resultant is fully determined,
    # i.e. the computation is exact or the residue modulo p^precision is
    # nonzero
    def determined?
      exact? || !to_i.zero?
    end

    # determined p-adic valuation of the resultant; Infinity when the
    # resultant is exactly zero and nil when the current precision cannot
    # determine it
    def valuation
      return padic_valuation(value, prime) if exact?

      determined? ? padic_valuation(to_i, prime) : nil
    end

    # certified lower bound for the p-adic valuation of the resultant
    def valuation_bound
      determined? ? valuation : precision
    end

    # number of p-adic digits to which the unit part of the resultant is
    # still known after factoring out the determined valuation; 0 when the
    # valuation itself is undetermined
    def unit_precision
      return Float::INFINITY if exact?

      determined? ? precision - valuation : 0
    end

    # true only when the resultant is known to be exactly zero
    def zero?
      exact? && value.zero?
    end

    # true if a common factor is certified, false if the resultant is
    # certified nonzero, and nil when the current precision cannot decide
    def common_factor?
      return true if zero?

      determined? ? false : nil
    end

    def to_s
      to_i.to_s
    end

    def inspect
      "<ResultantResult: #{self}#{precision_string}>"
    end

    private

    def precision_string
      exact? ? " (exact)" : " (mod #{prime}^#{precision})"
    end
  end
end
