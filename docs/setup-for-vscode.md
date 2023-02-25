# Set up your environment for vscode
* Install the "Modern Fortran" extension
* Install fortls
* Restart vscode
<br>
You can then debug using the "Debug CI_gen in example/test" launch configuration or pressing F5.

## fortls installation
Prerequisite: python3.7
```sh
sudo apt install python3.7
## Do this if you already have other python versions installed
sudo update-alternatives --install /usr/bin/python python /usr/bin/python3.7 1
sudo update-alternatives --config python # then choose 1
```
Note: you must install fortls using python3.7
```sh
python3.7 -m pip install fortls
```
fortls path should be: $HOME/.local/bin/fortls
Then update .vscode/settings.json with your path

## Additional notes
You may need to install this dependency:
```sh
sudo apt install liblapack-dev
```