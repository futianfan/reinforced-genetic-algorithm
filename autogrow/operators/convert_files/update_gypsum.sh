branch=$1
if [ "$branch" == "" ]; then
    echo "Branch will be whatever is default: "
fi


# Remove all previous gypsum_dl
rm -rf gypsum_dl/


git clone https://git.durrantlab.pitt.edu/jdurrant/gypsum_dl.git
cd gypsum_dl/
pwd
git status
echo ""
 
git pull
git fetch
git pull
echo "############"
git status
echo "############"
if [ "$branch" != "" ]; then
    git checkout $branch
fi

git pull
git fetch
git pull

echo "############"
git status
echo "############"


# Remove the .git files
rm -rf .git
rm -rf .gitignore
echo ""

# Remove the .git files from dimorphite
cd gypsum_dl/Steps/SMILES/dimorphite_dl/
ls -a
rm -rf .git
rm -rf .gitignore
echo ""

ls -a