# Third-Party Notices

This repository's workflow helper and installation code is licensed under the
MIT License in `LICENSE`.

The release package may include `archives/vendor.tar.gz`, an offline cache of
third-party source archives used only to make installation reproducible on
restricted networks. These third-party projects are not relicensed by this
repository. They remain governed by their upstream licenses and terms.

## Bundled Source Archives

| Component | Upstream | Pinned revision/version | License status |
|---|---|---:|---|
| SCENIC+ | https://github.com/aertslab/scenicplus | e82b82f14b76618b850dfe442efc2421bb34f3b4 | Academic non-commercial / proprietary license, SPDX `NOASSERTION`; see upstream `LICENCE.txt`. |
| pycisTopic | https://github.com/aertslab/pycisTopic | 219225df56b32738d82cd14532b187a1483de04f | Academic non-commercial license; see upstream `LICENSE.txt`. |
| pycistarget | https://github.com/aertslab/pycistarget | 5aa517604e4842539a7531c16905825dc7cb80fb | Academic non-commercial license; see upstream `LICENSE.txt`. |
| create_cisTarget_databases | https://github.com/aertslab/create_cisTarget_databases | 304d5dc1b15e5c923908a50a1ec291c3faaccf9c | No standard license file was present in the bundled archive or GitHub license metadata at packaging time. Verify upstream terms before redistribution beyond internal research use. |
| Cluster-Buster | https://github.com/ghuls/cluster-buster | 5911cd6201b767a43316ce613afc6c9255dc3511 | No standard license file was present in the bundled archive. Source comments indicate unrestricted use/redistribution intent by the original author; verify upstream terms before redistribution beyond internal research use. |
| hdWGCNA | https://github.com/smorabit/hdWGCNA | afa09abb890f5be087b63e510a7346e8e1952ecc | GPL-3.0 according to the bundled R `DESCRIPTION` file. |
| AutoZyme | https://github.com/ArcInstitute/autozyme | 35f91f2229eb44d82710470803865d3c15102716 | MIT according to bundled `LICENSE` and package metadata. |
| LoomXpy | https://github.com/aertslab/LoomXpy | 61995ff10940968eac2cee8fe48300ab477a15d0 | MIT according to bundled `LICENSE`. |
| MALLET | https://mimno.github.io/Mallet/ | 2.0.8 | Common Public License 1.0 or later according to bundled `LICENSE`. |

## Practical Use Notes

- The workflow helper and installation code in this repository is permissively licensed for reuse.
- The bundled SCENIC+/pycisTopic/pycistarget layer is not a general commercial
  open-source grant. It is intended here for academic/laboratory research
  installation workflows.
- If distributing modified release packages, keep this notice, the root
  `LICENSE`, and the original third-party source archives or their upstream
  license files together.
