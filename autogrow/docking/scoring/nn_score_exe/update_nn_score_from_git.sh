branch_NN1=$1
if [ "$branch_NN1" == "" ]; then
    echo "Branch will be whatever is default: "
fi
branch_NN2=$2
if [ "$branch_NN2" == "" ]; then
    echo "Branch will be whatever is default: "
fi


# Remove both NNScores
rm -rf NNScore_*/

git clone https://git.durrantlab.pitt.edu/jdurrant/nnscore1.git

# cp __init__.py nnscore1/

cd nnscore1/
pwd
git status
echo ""
 
git pull
git fetch
git pull
echo "############"
git status
echo "############"
if [ "$branch_NN1" != "" ]; then
    git checkout $branch_NN1
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
ls -a


cd ../

############################
#Update NN2
############################

git clone https://git.durrantlab.pitt.edu/jdurrant/nnscore2.git
# cp __init__.py nnscore2/

cd nnscore2/
pwd
git status
echo ""
 
git pull
git fetch
git pull
echo "############"
git status
echo "############"
if [ "$branch_NN2" != "" ]; then
    git checkout $branch_NN2
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
ls -a


cd ../