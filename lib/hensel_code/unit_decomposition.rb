# frozen_string_literal: true

module HenselCode
  class UnitDecomposition
    attr_reader :status, :sign_status, :prime, :exponent, :unit_exponent, :valuation,
                :valuation_power, :teichmuller_root, :principal_unit, :sign

    def initialize(owner, status, valuation = nil, teichmuller_root = nil, principal_unit = nil, sign = nil,
                   sign_status = nil)
      @owner = owner
      @status = status
      @prime = owner.prime
      @exponent = owner.exponent
      @valuation = valuation
      @valuation_power = valuation.nil? ? nil : prime**valuation
      @teichmuller_root = teichmuller_root
      @principal_unit = principal_unit
      @unit_exponent = teichmuller_root.nil? ? nil : teichmuller_root.exponent
      @sign = sign
      @sign_status = sign_status
    end

    def nonzero?
      status == :nonzero
    end

    alias state status
    alias omega teichmuller_root
    alias root teichmuller_root
    alias root_of_unity teichmuller_root
    alias root_unit teichmuller_root
    alias principal principal_unit
    alias valuation_component valuation_power

    def zero?
      status == :zero
    end

    def zero_mod_prime_power?
      status == :zero_mod_prime_power
    end

    alias known_zero_mod_prime_power? zero_mod_prime_power?

    def unit_modulus
      return nil if unit_exponent.nil?

      prime**unit_exponent
    end

    def root_power_equation_satisfied?
      return false unless nonzero?

      power = prime == 2 ? 2 : prime - 1
      teichmuller_root.residue.pow(power, unit_modulus) == 1
    end

    def principal_unit_condition_satisfied?
      return false unless nonzero?

      modulus = prime == 2 && unit_exponent >= 2 ? 4 : prime
      principal_unit.residue % modulus == 1
    end

    def compatible_with?(other)
      other = other.unit_decomposition if other.respond_to?(:unit_decomposition)
      return false unless other.is_a?(UnitDecomposition)
      return false unless prime == other.prime && status == other.status && valuation == other.valuation
      return status == :zero || status == :zero_mod_prime_power unless nonzero?

      precision = [unit_exponent, other.unit_exponent].min
      modulus = prime**precision
      teichmuller_root.residue % modulus == other.teichmuller_root.residue % modulus &&
        principal_unit.residue % modulus == other.principal_unit.residue % modulus
    end

    def recombine
      residue = if nonzero?
                  unit_residue = (teichmuller_root.residue * principal_unit.residue) % unit_modulus
                  (valuation_power * unit_residue) % (prime**exponent)
                else
                  0
                end
      @owner.send(:build_padic_integer, residue, exponent)
    end

    alias restore recombine

    def to_h
      {
        status: status,
        prime: prime,
        exponent: exponent,
        unit_exponent: unit_exponent,
        valuation: valuation,
        valuation_power: valuation_power,
        teichmuller_root: teichmuller_root,
        principal_unit: principal_unit,
        sign: sign,
        sign_status: sign_status
      }
    end

    alias to_hash to_h

    def [](key)
      to_h[key]
    end

    def inspect
      "<HenselCode::UnitDecomposition: #{status}, p=#{prime}, N=#{exponent}, " \
        "v=#{valuation.inspect}, sign=#{sign.inspect}>"
    end
  end
end
