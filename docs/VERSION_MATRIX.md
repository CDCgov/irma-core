# IRMA vs. IRMA-core Compatibility Matrix

This matrix can help with tracking which version of IRMA-core shipped with which version of [IRMA]. Using IRMA's pre-packaged [releases](https://github.com/CDCgov/irma/releases) is preferred.

| IRMA version | IRMA-core version | Notes                                                                                              |
| ------------ | ----------------- | -------------------------------------------------------------------------------------------------- |
| v1.3.1       | 0.6.2             | Better errors, zipped reader used in `preprocess`; `sampler` functionality exists but not exposed. |
| v1.3.0       | 0.5.1             | Unified preprocessing + mandatory IRMA-core; `trimmer` functionality exists but not exposed.       |

## Legacy tagging

| IRMA-core compatibility tag | IRMA-core version | Notes                                              |
| --------------------------- | ----------------- | -------------------------------------------------- |
| IRMA@v1.1.5                 | 0.1.4             |                                                    |
| IRMA@v1.1.4                 | 0.1.3             | Fix unstable dependencies for the latest nightly.  |
| IRMA@v1.1.2                 | 0.1.2             | Updated to maintain compatibility with `std::simd` |
| IRMA@v1.1.1                 | 0.1.1             |                                                    |
| IRMA@v1.1.0                 | 0.1.0             | Initial release.                                   |

[IRMA]: https://github.com/CDCgov/irma
