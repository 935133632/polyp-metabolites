# Throttled launcher for the extended analyses (max 8 python processes at a time). Run from the repository root.
$ErrorActionPreference = "Stop"
$root = "D:\research\Xingang Li\polyp_metabolites"; Set-Location $root
$max = 8
$jobs = @()
foreach ($k in "K63", "K31") { foreach ($m in "M1", "M2a", "M2b", "M2c") { $jobs += ,@("code/12_model_adjustment.py $k $m", "adj_${k}_$m") } }
$codes = "K57","K21","K29","K44","K59","K62","K64","K92","K80","K40","K22","K52","K76","K58","K20","K30","K56","K66","K25","K42","K83","K43","K26","CRC"
foreach ($c in $codes) { $jobs += ,@("code/13_specificity.py $c", "spec_$c") }
$running = @()
foreach ($j in $jobs) {
    while (($running | Where-Object { -not $_.HasExited }).Count -ge $max) { Start-Sleep -Seconds 20 }
    $running += Start-Process python -ArgumentList $j[0] -RedirectStandardOutput "results\$($j[1]).log" -RedirectStandardError "results\$($j[1]).err" -WindowStyle Hidden -PassThru
}
$running | Wait-Process
"queue finished" | Out-File "results\queue_done.txt"
