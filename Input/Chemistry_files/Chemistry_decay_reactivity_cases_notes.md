# Legacy Note: decay/reactivity chemistry cases

These files were created for the older concentration-decay chemistry workflow.
That workflow is no longer used on the active path.

At present, chemistry input files are not read by the active path.
`Vref` is read from the simulation input file instead.

So files such as:

- `Chemistry_decay_reactivity_case_1.txt`
- `Chemistry_decay_reactivity_case_2.txt`
- `Chemistry_decay_reactivity_case_3.txt`

no longer define different decay or reactivity strengths by themselves.

If you want different chemistry behavior now, it must come from the active
code-side mineral-volume law, not from `C0` or `k_decay` values in these input
files.
