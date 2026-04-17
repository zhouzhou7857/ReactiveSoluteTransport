# Legacy Note: reactivity chemistry cases

These files were created for the older concentration-decay chemistry workflow.
That workflow is no longer used on the active path.

At present, chemistry input files are not read by the active path.
`Vref` is read from the simulation input file instead.

So files such as:

- `Chemistry_reactivity_case_1.txt`
- `Chemistry_reactivity_case_2.txt`
- `Chemistry_reactivity_case_3.txt`

do not change reactivity strength unless the active code-side chemistry law is
also changed.
