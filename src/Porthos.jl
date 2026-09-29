"""
    Porthos

Port-Hamiltonian Operation and Stability: power system simulation, stability studies and
region-of-attraction certificates from one generic model code. See `README.md` and
`docs/ROADMAP.md`.
"""
module Porthos

using LinearAlgebra
using SparseArrays
using Printf: Printf
using ForwardDiff: ForwardDiff
import SciMLBase
import Sundials
import JLD2
using Dates: Dates
using Artifacts: Artifacts
using SHA: SHA
using Random: Random
import JSON3
import JSONSchema

# io: inputs, parameters, parity pack
include("io/expr.jl")
include("io/json.jl")
include("io/schema.jl")
include("io/case.jl")
include("io/scenario.jl")
include("io/contracts.jl")
include("io/params.jl")
include("io/parity.jl")

# network
include("network/ybus.jl")
include("network/events.jl")

# power flow
include("powerflow/newton.jl")

# components
include("components/primitives.jl")
include("components/interface.jl")
include("components/machines/genrou.jl")
include("components/machines/gensal.jl")
include("components/exciters/ieeet1.jl")
include("components/governors/ieeeg1.jl")
include("components/governors/ieeeg3.jl")
include("components/loads/complexload.jl")
include("components/dispatch.jl")
include("components/observables.jl")

# assembly
include("assembly/wiring.jl")
include("assembly/dae.jl")
include("assembly/sparsity.jl")

# initialisation
include("init/components.jl")
include("init/equilibrium.jl")

# simulation
include("sim/bdf1.jl")
include("sim/ida.jl")
include("sim/results.jl")

# port-Hamiltonian audits
include("ph/storage.jl")
include("ph/audit.jl")
include("ph/ports.jl")

# PowerFactory interface (runs pf/, reads its results)
include("io/powerfactory.jl")

# cached compilation of the simulation path
include("precompile.jl")

# io
export parse_param_expr, param_value
export load_case, write_case, load_scenario, write_scenario, load_json_input, validate_json
export Case, CaseConfig, BusData, LoadData, SourceData, LineData, ShuntData, ComponentSpec, Wire
export Scenario, SolverSettings, AbstractEvent, BusFault, LineFault, OtherEvent
export load_contracts, contract, contract_key, ContractSet, ContractEntry, DomainClause
export component, param, hasparam, component_params, MODEL_TYPES, check_model_type,
       UnsupportedModelError
export ParityPack, load_parity_pack, verify_parity_pack, parity_pack_dir
# network
export Network, nbus, bus_index, ybus, ybus_pf, ybus_dae, norton_stamps, NortonStamp
export load_admittances, LoadAdmittances, fault_admittance, fault_shunts, FaultShunt,
       with_fault, split_line_for_fault
# components
export AbstractComponent, build_component, with_params, model_type, state_names, input_names,
       output_names, nstates, ninputs, noutputs, params, param_dict, ports, bus, rhs!,
       outputs!, step_outputs!, modes, hamiltonian, grad_hamiltonian!, grad_hamiltonian,
       injection, norton_admittance, default_contracts
export GENROU_PHTRUE, GENSAL_PHTRUE, IEEET1_PHTRUE, IEEEG1_PHTRUE, IEEEG3_PHTRUE, COMPLEXLOAD
export NoModes, ModeLog, UndecidedBranch, component_role
# assembly
export DAESystem, assemble, dae_residual!, dae_residual, nalg, resolve_wiring, InputSource,
       jacobian_pattern
# initialisation
export init_from_phasor, init_from_targets, MachineTargets, first_pass, solve_equilibrium,
       EquilibriumResult, component_io
# simulation
export DAEWorkspace, SimResult, simulate_bdf1, simulate_ida, consistent_voltages!,
       simulate, csv_columns, csv_row, write_results_csv, write_results_jld2, run_metadata,
       observable_names, observable_values!
# port-Hamiltonian audits
export storage_components, total_hamiltonian, grad_total_hamiltonian, hessian_total_hamiltonian,
       solve_network, reduced_field, reduced_jacobian, PhysicalProjection, reservoir_states,
       physical_projection, shifted_storage_audit, PortModel, port_model, transfer,
       real_part_crossings, passivity_certificate, port_zeros, frequency_response,
       loop_port_model
# PowerFactory
export PFResults, read_pf_results, pf_signal, pf_command, pf_simulate, pf_inspect,
       pf_machine_map, pf_compare
# power flow
export solve_powerflow, PowerFlowSpec, PowerFlowResult, bus_power, BusType, PQ_BUS, PV_BUS,
       SLACK_BUS

end
