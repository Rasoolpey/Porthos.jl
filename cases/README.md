# Cases

Case (system) and scenario JSON in the PHPS format, copied unchanged from PHPS_Opt
`phps/cases/` at commit `ba11ea1827c344ace9dd2734215b5f5b4dbbe822` (the committed blobs,
byte for byte). Only the directories the parity cases need are here:

| Directory | Parity case |
|---|---|
| `IEEE39Bus_PF/` | `base` (`system_phtrue.json`, `bus_fault_bus16_150ms.json`) |
| `IEEE39Bus_PF_gfl/` | `gfl` |
| `IEEE39Bus_PF_gfl-zif/` | `gfl_zif` |
| `IEEE39Bus_PF_gfm-vsm/` | `vsm` |
| `IEEE39Bus_PF_gfm-droop/` | `droop` |
| `IEEE39Bus_PF_gfm-voc/` | `voc` |

`schema/` holds the JSON Schemas (`system.schema.json`, `scenario.schema.json`) that
`load_case` and `load_scenario` validate against. They accept every one of the 200 JSON
files in PHPS `phps/cases` (checked by the test suite when a PHPS checkout is present).

Do not edit these files by hand. A case variant is a new file (JSON edits are what the
Python layer is allowed to do); a newer PHPS copy is recorded here with its commit.
