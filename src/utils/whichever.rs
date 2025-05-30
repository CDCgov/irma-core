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
            $(
                $(#[$variant_meta:meta])*
                $variant:ident($ty:ty)
            ),+
            $(,)?
        }
        $(impl $trait:ident for $struct_name2:ident {$($impl_block:tt)*}),*
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
                $(#[$variant_meta])*
                $variant($ty),
            )+
        }

        $(define_whichever!(@impl_trait $struct_name2, $trait, $($impl_block)*);)*
    };

    (@impl_trait $struct_name:ident, Read, $($impl_block:tt)*) => {
        impl ::std::io::Read for $struct_name {
            $($impl_block)*

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

    (@impl_trait $struct_name:ident, Write, $($impl_block:tt)*) => {
        impl ::std::io::Write for $struct_name {
            $($impl_block)*

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

    (@impl_trait $struct_name:ident, Iterator, $($impl_block:tt)*) => {
        impl Iterator for $struct_name {
            $($impl_block)*

            #[inline]
            fn next(&mut self) -> Option<Self::Item> {
                match_macro!(self, inner => inner.next())
            }

            #[inline]
            fn size_hint(&self) -> (usize, Option<usize>) {
                match_macro!(self, inner => inner.size_hint())
            }

            #[inline]
            fn count(self) -> usize {
                match_macro!(self, inner => inner.count())
            }

            #[inline]
            fn last(self) -> Option<Self::Item> {
                match_macro!(self, inner => inner.last())
            }

            #[inline]
            fn nth(&mut self, n: usize) -> Option<Self::Item> {
                match_macro!(self, inner => inner.nth(n))
            }

            #[inline]
            fn for_each<F>(self, f: F)
            where
                Self: Sized,
                F: FnMut(Self::Item),
            {
                match_macro!(self, inner => inner.for_each(f))
            }

            #[inline]
            fn collect<B>(self) -> B
            where
                B: FromIterator<Self::Item>,
                Self: Sized,
            {
                match_macro!(self, inner => inner.collect())
            }

            #[inline]
            fn partition<B, F>(self, f: F) -> (B, B)
            where
                Self: Sized,
                B: Default + Extend<Self::Item>,
                F: FnMut(&Self::Item) -> bool,
            {
                match_macro!(self, inner => inner.partition(f))
            }

            #[inline]
            fn try_fold<B, F, R>(&mut self, init: B, f: F) -> R
            where
                Self: Sized,
                F: FnMut(B, Self::Item) -> R,
                R: ::std::ops::Try<Output = B>,
            {
                match_macro!(self, inner => inner.try_fold(init, f))
            }

            #[inline]
            fn try_for_each<F, R>(&mut self, f: F) -> R
            where
                Self: Sized,
                F: FnMut(Self::Item) -> R,
                R: ::std::ops::Try<Output = ()>,
            {
                match_macro!(self, inner => inner.try_for_each(f))
            }

            #[inline]
            fn fold<B, F>(self, init: B, f: F) -> B
            where
                Self: Sized,
                F: FnMut(B, Self::Item) -> B,
            {
                match_macro!(self, inner => inner.fold(init, f))
            }

            #[inline]
            fn reduce<F>(self, f: F) -> Option<Self::Item>
            where
                Self: Sized,
                F: FnMut(Self::Item, Self::Item) -> Self::Item,
            {
                match_macro!(self, inner => inner.reduce(f))
            }

            #[inline]
            fn all<F>(&mut self, f: F) -> bool
            where
                Self: Sized,
                F: FnMut(Self::Item) -> bool,
            {
                match_macro!(self, inner => inner.all(f))
            }

            #[inline]
            fn any<F>(&mut self, f: F) -> bool
            where
                Self: Sized,
                F: FnMut(Self::Item) -> bool,
            {
                match_macro!(self, inner => inner.any(f))
            }

            #[inline]
            fn find<P>(&mut self, predicate: P) -> Option<Self::Item>
            where
                Self: Sized,
                P: FnMut(&Self::Item) -> bool,
            {
                match_macro!(self, inner => inner.find(predicate))
            }

            #[inline]
            fn find_map<B, F>(&mut self, f: F) -> Option<B>
            where
                Self: Sized,
                F: FnMut(Self::Item) -> Option<B>,
            {
                match_macro!(self, inner => inner.find_map(f))
            }

            #[inline]
            fn position<P>(&mut self, predicate: P) -> Option<usize>
            where
                Self: Sized,
                P: FnMut(Self::Item) -> bool,
            {
                match_macro!(self, inner => inner.position(predicate))
            }

            #[inline]
            fn max_by_key<B, F>(self, f: F) -> Option<Self::Item>
            where
                B: Ord,
                Self: Sized,
                F: FnMut(&Self::Item) -> B,
            {
                match_macro!(self, inner => inner.max_by_key(f))
            }

            #[inline]
            fn max_by<F>(self, compare: F) -> Option<Self::Item>
            where
                Self: Sized,
                F: FnMut(&Self::Item, &Self::Item) -> ::std::cmp::Ordering,
            {
                match_macro!(self, inner => inner.max_by(compare))
            }

            #[inline]
            fn min_by_key<B, F>(self, f: F) -> Option<Self::Item>
            where
                B: Ord,
                Self: Sized,
                F: FnMut(&Self::Item) -> B,
            {
                match_macro!(self, inner => inner.min_by_key(f))
            }

            #[inline]
            fn min_by<F>(self, compare: F) -> Option<Self::Item>
            where
                Self: Sized,
                F: FnMut(&Self::Item, &Self::Item) -> ::std::cmp::Ordering,
            {
                match_macro!(self, inner => inner.min_by(compare))
            }

            #[inline]
            fn sum<S>(self) -> S
            where
                Self: Sized,
                S: ::std::iter::Sum<Self::Item>,
            {
                match_macro!(self, inner => inner.sum())
            }

            #[inline]
            fn product<S>(self) -> S
            where
                Self: Sized,
                S: ::std::iter::Product<Self::Item>,
            {
                match_macro!(self, inner => inner.product())
            }

            #[inline]
            fn is_sorted_by<F>(self, compare: F) -> bool
            where
                Self: Sized,
                F: FnMut(&Self::Item, &Self::Item) -> bool,
            {
                match_macro!(self, inner => inner.is_sorted_by(compare))
            }

            #[inline]
            fn is_sorted_by_key<F, K>(self, f: F) -> bool
            where
                Self: Sized,
                F: FnMut(Self::Item) -> K,
                K: PartialOrd,
            {
                match_macro!(self, inner => inner.is_sorted_by_key(f))
            }
        }
    };

    (@impl_trait $struct_name:ident, $other:ident) => {
        compile_error!(concat!("Unsupported trait: ", stringify!($other)));
    };
}

pub(crate) use define_whichever;
