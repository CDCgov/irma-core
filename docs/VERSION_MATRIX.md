# IRMA vs. IRMA-core Compatibility Matrix

This matrix can help with tracking which version of IRMA-core shipped with which version of [IRMA]. Using IRMA's pre-packaged [releases](https://github.com/CDCgov/irma/releases) is preferred.

| IRMA version              | IRMA-core version | Notes                                                                                           |
| ------------------------- | ----------------- | ----------------------------------------------------------------------------------------------- |
| [0.5.0-irma-compat-1.3.0] | 0.5.0             | Unified preprocessing + mandatory IRMA-core; new trimming functionality exists but not exposed. |
| [IRMA@v1.1.5]             | 0.1.4             |                                                                                                 |
| [IRMA@v1.1.4]             | 0.1.3             | Fix unstable dependencies for the latest nightly.                                               |
| [IRMA@v1.1.2]             | 0.1.2             | Updated to maintain compatibility with `std::simd`                                              |
| [IRMA@v1.1.1]             | 0.1.1             |                                                                                                 |
| [IRMA@v1.1.0]             | 0.1.0             | Initial release.                                                                                |

[0.5.0-irma-compat-1.3.0]: https://github.com/CDCgov/irma-core/releases/tag/0.5.0-irma-compat-1.3.0
[IRMA@v1.1.5]: https://github.com/CDCgov/irma-core/releases/tag/IRMA@v1.1.5
[IRMA@v1.1.4]: https://github.com/CDCgov/irma-core/releases/tag/IRMA@v1.1.4
[IRMA@v1.1.2]: https://github.com/CDCgov/irma-core/releases/tag/IRMA@v1.1.2
[IRMA@v1.1.1]: https://github.com/CDCgov/irma-core/releases/tag/IRMA@v1.1.1
[IRMA@v1.1.0]: https://github.com/CDCgov/irma-core/releases/tag/IRMA@v1.1.0
[IRMA]: https://github.com/CDCgov/irma
