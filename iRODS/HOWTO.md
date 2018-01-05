# How to configure your local environment for iRODS

1. Create an `~/irods` folder and copy
   [irods_environment.json](irods_environment.json) inside it. This file
   specifies the server settings.

2. _If not already installed_, clone [ishell](https://github.com/niess/ishell)
    and do a local install by following the installation instructions. Then,
    in order to add the ishell binaries to your path you can source
    the `ishell/setup.sh` file.

    On **fh1** ishell is _already installed_ in `/project/fh1-project-huepra/le6232/soft/ishell`. Though you'll need to
    source the `setup.sh` file or add this step to your .bashrc.

3. Run `iinit`, provided by ishell. Note that you need to know our iRODS
   password.

This procedure needs to be done only once. Once properly configured you can
navigate the iRODS data by issuing the `ìshell` command. Alternatively you can
also execute ishell scripts, see for example [put-events.ish](put-events.ish),
or run ishell in interpreted mode with the `-c` option, e.g.
`ishell -c "cd grand; ls"`.
