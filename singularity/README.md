# singularity

Definition files for specific arch for applications and 
general apporaches ("sra", "alignment", etc).

## Sherlock and building images

Building singularity images requires `root` so this cannot be done
on `sherlock`.
My local fix on M1 mac is emulating `centOS8` and building there.
The resulting `.sif` files are too large to track with GH in a 
scalable manner.


