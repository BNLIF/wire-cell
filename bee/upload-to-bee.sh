LOGIN_URL=https://www.phy.bnl.gov/twister/bee
COOKIES=cookies.txt
CURL_BIN="curl -k -s -c $COOKIES -b $COOKIES -e $LOGIN_URL"

$CURL_BIN $LOGIN_URL > /dev/null
echo -n "Django Auth: get csrftoken ... "
DJANGO_TOKEN="$(grep csrftoken $COOKIES | sed 's/^.*csrftoken\s*//')"
echo $DJANGO_TOKEN

FILENAME=$1
echo "uploading file $FILENAME ... "
UUID="$($CURL_BIN \
    -H "X-CSRFToken: $DJANGO_TOKEN" \
    -F "file=@$FILENAME" \
    $LOGIN_URL/upload/)"

# echo "$DJANGO_TOKEN&file=@/home/vagrant/Sites/bee/test/test.zip"
# echo "Please redirect your browser to:"
URL="$LOGIN_URL/set/$UUID/event/list/"
${BROWSER:-echo} "$URL"

rm $COOKIES

