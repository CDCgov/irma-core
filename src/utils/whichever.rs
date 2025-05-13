/// A macro to define an enum similar to Either, but with any number of
/// variants, each containing a known (not generic) type. An invocation will
/// contain:
/// 1. The enum definition, with any number of outer attributes, an optional
///    visibility specifier, and the variants
/// 2. The traits to implement, using `impl Trait {}`. Currently, we support
///    `Read` and `Write`. The braces should be left empty (the implementations
///    will be filled automatically).
macro_rules! define_whichever {
    (
        $(#[$meta:meta])*
        $vis:vis enum $struct_name:ident {
            $($variant:ident($ty:ty)),+
            $(,)?
        }
        $(impl $trait:ident for $struct_name2:ident {}),*
    ) => {
        macro_rules! match_macro {
            ($value:expr, $pattern:pat => $result:expr) => {
                match $value {
                    $(
                        $struct_name::$variant($pattern) => $result,
                    )+
                }
            };
        }

        $(#[$meta])*
        $vis enum $struct_name {
            $(
                $variant($ty),
            )+
        }

        $(define_whichever!(@impl_trait $struct_name2, $trait);)*
    };

    (@impl_trait $struct_name:ident, Read) => {
        impl ::std::io::Read for $struct_name {
            #[inline]
            fn read(&mut self, buf: &mut [u8]) -> ::std::io::Result<usize> {
                match_macro!(self, inner => inner.read(buf))
            }

            #[inline]
            fn read_vectored(&mut self, bufs: &mut [::std::io::IoSliceMut<'_>]) -> ::std::io::Result<usize> {
                match_macro!(self, inner => inner.read_vectored(bufs))
            }

            #[inline]
            fn read_to_end(&mut self, buf: &mut ::std::vec::Vec<u8>) -> ::std::io::Result<usize> {
                match_macro!(self, inner => inner.read_to_end(buf))
            }

            #[inline]
            fn read_to_string(&mut self, buf: &mut ::std::string::String) -> ::std::io::Result<usize> {
                match_macro!(self, inner => inner.read_to_string(buf))
            }

            #[inline]
            fn read_exact(&mut self, buf: &mut [u8]) -> ::std::io::Result<()> {
                match_macro!(self, inner => inner.read_exact(buf))
            }
        }
    };

    (@impl_trait $struct_name:ident, Write) => {
        impl ::std::io::Write for $struct_name {
            #[inline]
            fn write(&mut self, buf: &[u8]) -> ::std::io::Result<usize> {
                match_macro!(self, inner => inner.write(buf))
            }

            #[inline]
            fn flush(&mut self) -> ::std::io::Result<()> {
                match_macro!(self, inner => inner.flush())
            }

            fn write_vectored(&mut self, bufs: &[::std::io::IoSlice<'_>]) -> ::std::io::Result<usize> {
                match_macro!(self, inner => inner.write_vectored(bufs))
            }

            fn write_all(&mut self, buf: &[u8]) -> ::std::io::Result<()> {
                match_macro!(self, inner => inner.write_all(buf))
            }

            fn write_fmt(&mut self, fmt: ::std::fmt::Arguments<'_>) -> ::std::io::Result<()> {
                match_macro!(self, inner => inner.write_fmt(fmt))
            }
        }
    };

    (@impl_trait $struct_name:ident, $other:ident) => {
        compile_error!(concat!("Unsupported trait: ", stringify!($other)));
    };
}

pub(crate) use define_whichever;
