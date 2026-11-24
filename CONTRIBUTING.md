# Contributing

Thank you for helping. Problems and questions: open an issue on GitHub, with the command you ran, the full error
message and the output of `pixi --version`.

To change code:

1. Make a branch: `git switch -c short-description`.
2. Set up and check: `pixi install --locked`, then `pixi run test` (pytest, nf-test, stub run) and
   `pixi run -e lint lint`. All must pass.
3. Add software only with `pixi add` (commit `pixi.toml` and `pixi.lock` together); put settings in
   `config/params.yaml`, never paths of your machine (those go in `config/site.config`, which is not committed).
4. Commit in small steps with a type prefix: `add:`, `fix:`, `refactor:`, `docs:`, for example
   `fix: keep salmon's NumReads as floats`.
5. Open a pull request. The *Tests* and *Reproduce* checks must be green before it is merged.

Everyone taking part follows the [Code of Conduct](CODE_OF_CONDUCT.md).
