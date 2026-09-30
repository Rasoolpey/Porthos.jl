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
import IntervalArithmetic

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
include("components/converters/common.jl")
include("components/converters/gfm_vsm.jl")
include("components/converters/gfm_droop.jl")
include("components/converters/gfm_voc.jl")
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
include("ph/dissipativity.jl")
include("ph/terminal.jl")
include("ph/structure.jl")
include("ph/power.jl")

# ROA certificates
include("roa/interval.jl")
include("roa/section.jl")
include("roa/candidate.jl")
include("roa/enclosure.jl")
include("roa/containment.jl")
include("roa/certificate.jl")
include("roa/analytic.jl")
include("roa/check.jl")

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
export GENROU_PHTRUE, GENSAL_PHTRUE, IEEET1_PHTRUE, IEEEG1_PHTRUE, IEEEG3_PHTRUE, COMPLEXLOAD,
       GFM_VSM_PHTRUE, GFM_DROOP_PHTRUE, GFM_VOC_PHTRUE
export NoModes, ModeLog, UndecidedBranch, component_role
# assembly
export DAESystem, assemble, dae_residual!, dae_residual, nalg, resolve_wiring, InputSource,
       jacobian_pattern
# initialisation
export init_from_phasor, init_from_targets, MachineTargets, first_pass, solve_equilibrium,
       converter_init, converter_current, lag_states,
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
       loop_port_model, MultiPortModel, open_loops_model, multiport_margin, kyp_riccati,
       port_storage, rest_storage, port_margin, loop_margin, TerminalModel, NetworkModel,
       sync_jacobian, terminal_models, terminal_margins, state_groups, storage_pattern,
       reference_section, lyapunov_check, structured_lyapunov, pow2_scaling, decay_margin,
       margin_residuals, verified_min_eig, pattern_certificate, verified_lyapunov,
       rank_couplings, add_coupling!, couple_states!, section_pattern,
       component_power, network_power, power_audit_samples, port_power_audit
# ROA certificates
export ProofFailure, is_interval_type, interval_solve, verified_max_eig, weyl_max_eig, verified_inverse_diagonal,
       ellipsoid_half_widths, SectionModel, section_model, lift, section_coordinates,
       section_field, state_field, project_field, section_residual, section_jacobian,
       LyapunovCandidate, QuadraticCandidate, quadratic_candidate, candidate_model,
       candidate_value, candidate_gradient, positivity_proof, sublevel_half_widths,
       gradient_matrix_hull, candidate_fingerprint, candidate_record, model_fingerprint,
       EquilibriumEnclosure, enclose_equilibrium, KCLBranch, enclose_kcl_branch, krawczyk,
       jacobian_hull, centered_hull, ContainmentAudit, containment_audit, box_digest, certify_level,
       certify_roa, section_record, software_record, write_certificate, eval_param_expr,
       system_digest, interval_cholesky, cholesky_positive_definite, roa_check,
       AbstractSectionModel, AnalyticModel, rotation_action, check_rotation_symmetry,
       certificate_claim, model_assumptions
# PowerFactory
export PFResults, read_pf_results, pf_signal, pf_command, pf_simulate, pf_inspect,
       pf_machine_map, pf_compare
# power flow
export solve_powerflow, PowerFlowSpec, PowerFlowResult, bus_power, BusType, PQ_BUS, PV_BUS,
       SLACK_BUS

end
