#ifndef CELL_HPP
#define CELL_HPP

/**
 * @class Cellmodel
 * @brief Abstract base class representing a generic cell model for simulations.
 */
class Cellmodel {
   protected:
    /**
     * @brief Protected constructor to prevent direct instantiation.
     */
    Cellmodel() {}

   public:
    unsigned short algebraic_size;        ///< Size of the algebraic array
    unsigned short constants_size;        ///< Size of the constants array
    unsigned short states_size;           ///< Size of the states array
    unsigned short gates_size;            ///< Size of the gates array
    unsigned short current_size;          ///< Size of the current array
    unsigned short concs_size;            ///< Size of the concentrations array
    double ALGEBRAIC[255];                ///< Array to store algebraic values
    double CONSTANTS[255];                ///< Array to store constants
    double RATES[255];                    ///< Array to store rates
    double STATES[255];                   ///< Array to store states
    char gates_header[255];               ///< Header for gates
    unsigned short gates_indices[255];    ///< Indices for gates
    char current_header[255];             ///< Header for current
    unsigned short current_indices[255];  ///< Indices for current
    char concs_header[255];               ///< Header for concentrations
    unsigned short concs_indices[255];    ///< Indices for concentrations

    /**
     * @brief Virtual destructor for proper cleanup of derived classes.
     */
    virtual ~Cellmodel() {}

    /**
     * @brief Pure virtual function to initialize constants.
     */
    virtual void initConsts() = 0;

    /**
     * @brief Initialize constants with a given type.
     * @param type Type of the cell model
     */
    virtual void initConsts(double type) {}

    /**
     * @brief Initialize constants with a given type, concentration, and Hill coefficients.
     * @param type Type of the cell model
     * @param conc Concentration value
     * @param hill Array of Hill coefficients
     */
    virtual void initConsts(double type, double conc, double *hill) {}

    /**
     * @brief Initialize constants with a given type, concentration, Hill coefficients, and Dutta flag.
     * @param type Type of the cell model
     * @param conc Concentration value
     * @param hill Array of Hill coefficients
     * @param is_dutta Flag indicating if the Dutta model is used
     */
    virtual void initConsts(double type, double conc, double *hill, bool is_dutta) {}

    /**
     * @brief Pure virtual function to compute rates based on time, constants, rates, states, and algebraic values.
     * @param TIME Time value
     * @param CONSTANTS Array of constants
     * @param RATES Array of rates
     * @param STATES Array of states
     * @param ALGEBRAIC Array of algebraic values
     */
    virtual void computeRates(double TIME, double *CONSTANTS, double *RATES, double *STATES, double *ALGEBRAIC) = 0;

    /**
     * @brief Solve the model analytically over a given time step.
     * @param dt Time step value
     */
    virtual void solveAnalytical(double dt) {}
};

#endif  // CELL_HPP
