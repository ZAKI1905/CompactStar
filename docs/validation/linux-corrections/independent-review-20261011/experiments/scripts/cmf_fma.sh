cd /mnt/sdd/zaki/CompactStar/linux-transition-20261010
echo "=== CMF per-star rows (Linux Debug) ==="
grep -E "structure unchanged|B within" toolchain/logs/CompactStar-Debug-suite-remainder.log | sed 's/^[[:space:]]*//' | cut -c1-150
echo "=== (Linux Release) ==="
grep -E "structure unchanged|B within" toolchain/logs/CompactStar-Release-suite-remainder.log | sed 's/^[[:space:]]*//' | cut -c1-150
echo
echo "=== FMA audit of Linux CompactStar objects ==="
for m in Debug Release; do
  cc=toolchain/build/CompactStar-$m/compile_commands.json
  n=$(grep -c '"command"' $cc 2>/dev/null); off=$(grep -c 'ffp-contract=off' $cc 2>/dev/null); ma=$(grep -c 'march=' $cc 2>/dev/null)
  lib=$(find toolchain/build/CompactStar-$m -maxdepth 1 -name 'libCompactStar.a' | head -1)
  fma=$(objdump -d $lib 2>/dev/null | grep -c -E 'vfmadd|vfmsub|vfnmadd|vfnmsub')
  echo "$m: TUs=$n contract_off=$off march=$ma lib=$lib fma_instructions=$fma"
done
