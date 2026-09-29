"""
    Porthos

Port-Hamiltonian Operation and Stability: power system simulation, stability studies and
region-of-attraction certificates from one generic model code. See `README.md` and
`docs/ROADMAP.md`.
"""
module Porthos

using LinearAlgebra
using SparseArrays
using Artifacts: Artifacts
using SHA: SHA
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

# io
export parse_param_expr, param_value
export load_case, write_case, load_scenario, write_scenario, load_json_input, validate_json
export Case, CaseConfig, BusData, LoadData, SourceData, LineData, ShuntData, ComponentSpec, Wire
export Scenario, SolverSettings, AbstractEvent, BusFault, LineFault, OtherEvent
export load_contracts, contract, contract_key, ContractSet, ContractEntry, DomainClause
export component, param, hasparam, component_params
export ParityPack, load_parity_pack, verify_parity_pack, parity_pack_dir
# network
export Network, nbus, bus_index, ybus, ybus_pf, ybus_dae, norton_stamps, NortonStamp
export load_admittances, LoadAdmittances, fault_admittance, fault_shunts, FaultShunt,
       with_fault, split_line_for_fault
# power flow
export solve_powerflow, PowerFlowSpec, PowerFlowResult, bus_power, BusType, PQ_BUS, PV_BUS,
       SLACK_BUS

end
